
import logging
from time import perf_counter, gmtime

from heatsource9.setup.model_setup import ModelSetup
from heatsource9.model.model_routine import ModelRoutine
from heatsource9.io.console import print_console
from heatsource9.io.control_file import validate_control_file
from heatsource9.domain.clock import pretty_time, time_parts
from heatsource9.io.output_files import OutputWriter

logger = logging.getLogger(__name__)


class ModelRunner:
    """
    Coordinates building and running of a model simulation including
    model setup (via ModelSetup), time keeping (via Clock), 
    running the model routine (via ModelRoutine), and writing 
    outputs (via OutputWriter).
    """

    def __init__(self, control_file_path):
        self.control_file_path = validate_control_file(control_file_path)

        self.output = OutputWriter()
        self.model_setup = None
        self.model_routine = None

    def temperature(self):
        return self._run("temperature")

    def solar(self):
        return self._run("solar")

    def hydraulics(self):
        return self._run("hydraulics")

    def _run(self, run_type):
        """Build and execute a model run.

        Builds the simulation, binds the output writer to model object, and loops
        through each timestep in the simulation clock. At each step it advances
        model calculations, and writes outputs.
        """
        try:
            self.model_setup = ModelSetup(control_file_path=self.control_file_path, run_type=run_type)
            self.model_routine = ModelRoutine(self.model_setup.run_type)
            msg = "Building model"
            print_console(msg)
            msg = f"Heat Source Version: {self.model_setup.params['version']}"
            print_console(msg)
            simulation = self.model_setup.build()

            self.output.bind(
                reach=self.model_setup.reach,
                params=self.model_setup.params,
                run_type=self.model_setup.run_type,
            )

            msg = "Starting model run"
            print_console(msg)
            start = perf_counter()
            total_steps = self._calculate_step_count()

            # Track cumulative downstream outflow (cms) summed by timestep.
            outflow_total = 0.0
            if total_steps > 0:
                print_console("Timesteps", progress_bar=True, current=0, total=total_steps)

            # Initialize daily solar accumulators before first step.
            self._reset_daily_solar_accumulators(simulation, 0, 0, 0, force=True)
            for step_number, epoch_seconds in enumerate(simulation.clock):
                hour, minute, second, doy, jc = time_parts(epoch_seconds)
                self._reset_daily_solar_accumulators(simulation, hour, minute, second)
                self._step(simulation, epoch_seconds, hour, minute, second, doy, jc)

                # Add the current discharge (cms) at the mouth to the outflow total.
                outflow_total += self._mouth_flow(simulation)
                if total_steps > 0:
                    print_console("Timesteps", progress_bar=True, current=step_number + 1, total=total_steps)

            elapsed = perf_counter() - start
            msg = "Model run complete in %.2f seconds" % elapsed
            print_console(msg)

            # Sum node inflow (cms) across the simulation.
            inflow_total = self._total_inflow(simulation)
            # Report water balance as total inflow / total outflow (cms).
            msg = f"Water Balance (cms): {inflow_total:0.4f}/{outflow_total:0.4f}"
            print_console(msg)
            msg = f"Outputs: {self.model_setup.params['outputdir']}"
            print_console(msg)
        except Exception:
            msg = f"Error in {run_type.title()} run."
            logger.exception(msg)
            raise

    def _step(self, simulation, time, hour, minute, second, doy, jc):
        try:
            self.model_routine.advance_timestep(simulation.nodes, time, hour, minute, second, doy, jc)
            if self._outputdt_due(time):
                self.output.queue_dt_outputs(time)
            if self._day_complete(time):
                if self.model_setup.run_type in ("temperature", "solar"):
                    self.output.queue_daily_outputs(time)
                self.output.write_to_csv()
            if self._run_complete(time):
                self.output.write_to_csv()
                self.output.close()
        except Exception:
            logger.exception(
                "Model failed at model date/time=%s.",
                pretty_time(time),
            )
            raise

    def _outputdt_due(self, time):
        """Return True when the current timestep should be added to the non daily output write queue."""
        params = self.model_setup.params
        if time < params["modelstart"] or time >= params["modelend"]:
            return False
        return (int(time) - int(params["modelstart"])) % int(params["outputdt"]) == 0

    def _day_complete(self, time):
        """Return True when time is the final modeled timestep of the current day."""
        params = self.model_setup.params
        if time < params["modelstart"] or time >= params["modelend"]:
            return False
        current = gmtime(float(time))
        nxt = gmtime(float(time + params["dt"]))
        return (current.tm_year, current.tm_yday) != (nxt.tm_year, nxt.tm_yday)

    def _run_complete(self, time):
        """Return True when the current timestep is at or past modelend time."""
        return time >= float(self.model_setup.params["modelend"])

    def _reset_daily_solar_accumulators(
        self,
        simulation,
        current_hour,
        current_minute,
        current_second,
        force=False):
        """
        Reset F_DailySum and Solar_Blocked values. They are reset daily.
        """
        params = self.model_setup.params
        trans_count = int(params.get("trans_count", 0) or 0)
        transsample_count = int(params.get("transsample_count", 0) or 0)

        if not force:
            if not (
                current_hour == 0 and
                current_minute == 0 and
                current_second == 0
            ):
                return

        for nd in simulation.nodes:
            nd.F_DailySum = [0] * 5
            nd.Solar_Blocked = {i: [0] * transsample_count for i in range(trans_count)}
            nd.Solar_Blocked["diffuse"] = 0

    def _calculate_step_count(self):
        """Calculate number of time step iterations for progress reporting."""
        params = self.model_setup.params
        try:
            start = float(params["flushtimestart"])
            end = float(params["modelend"])
            dt = float(params["dt"])
        except (KeyError, TypeError, ValueError):
            return 0
        if dt <= 0 or end < start:
            return 0
        return int((end - start) // dt) + 1

    @staticmethod
    def _mouth_flow(simulation):
        """Return current discharge at the downstream node in cms."""
        mouth = simulation.nodes[-1]
        return float(mouth.Q)

    @staticmethod
    def _total_inflow(simulation):
        """Return total accumulated node inflow in cms."""
        return float(sum(float(node.Q_mass) for node in simulation.nodes))
