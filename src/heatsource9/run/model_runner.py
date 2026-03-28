"""
Model Run coordi

This mcoordinates:
  - ModelSetup (setup)
  - Clock (domain)
  - ModelRoutine (model)
  - OutputWriter (io)
"""

import logging
from pathlib import Path
from time import perf_counter

from heatsource9.setup.model_setup import ModelSetup
from heatsource9.model.model_routine import ModelRoutine
from heatsource9.io.console import print_console
from heatsource9.io.control_file import validate_control_file
from heatsource9.domain.clock import pretty_time, time_parts
from heatsource9.io.output_files import OutputWriter
from heatsource9.io.logging_config import configure_logging

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
        model_dir = Path(self.control_file_path).parent
        try:
            configure_logging(model_dir=model_dir, overwrite=True)
        except PermissionError as exc:
            raise PermissionError(
                f"Cannot write log file '{model_dir / 'heatsource.log'}'. Check write permissions for this folder."
            ) from exc

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
        self.model_setup = ModelSetup(control_file_path=self.control_file_path, run_type=run_type)
        self.model_routine = ModelRoutine(self.model_setup.run_type)
        logger.info("Building model")
        simulation = self.model_setup.build()

        self.output.bind(
            reach=self.model_setup.reach,
            params=self.model_setup.params,
            run_type=self.model_setup.run_type,
        )

        logger.info("Starting model run")
        start = perf_counter()
        total_steps = self._calculate_step_count()

        # Track cumulative downstream outflow (cms) summed by timestep.
        outflow_total = 0.0
        if total_steps > 0:
            print_console("Timesteps", progress_bar=True, current=0, total=total_steps)

        # Initialize daily solar accumulators before first step.
        self._reset_daily_solar_accumulators(simulation, 0, 0, 0, force=True)
        for step_number, epoch_seconds in enumerate(simulation.clock):
            hour, minute, second, jd, jc = time_parts(epoch_seconds)
            self._reset_daily_solar_accumulators(simulation, hour, minute, second)
            logger.debug("Step %d, %s", step_number, epoch_seconds)
            self._step(simulation, epoch_seconds, hour, minute, second, jd, jc)

            # Add the current discharge (cms) at the mouth to the outflow total.
            outflow_total += self._mouth_flow(simulation)
            if total_steps > 0:
                print_console("Timesteps", progress_bar=True, current=step_number + 1, total=total_steps)

        elapsed = perf_counter() - start
        logger.info("Model run complete in %.2f seconds", elapsed)

        # Sum node inflow (cms) across the simulation.
        inflow_total = self._total_inflow(simulation)

        # Report water balance as total inflow / total outflow (cms).
        print_console(f"Water Balance (cms): {inflow_total:0.4f}/{outflow_total:0.4f}")

    def _step(self, simulation, time, hour, minute, second, jd, jc):
        try:
            self.model_routine.advance_timestep(simulation.nodes, time, hour, minute, second, jd, jc)
            if minute == 0 and second == 0:
                self.output.write_step(simulation, time, hour, minute, second)
        except Exception:
            logger.exception(
                "Model step failed at model date/time=%s.",
                pretty_time(time),
            )
            raise

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
