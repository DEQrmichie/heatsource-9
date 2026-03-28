

class ModelRoutine:
    """
    Class that stores the run type and advances all stream nodes by one timestep 
    using the appropriate model run calculations.
    """

    def __init__(self, run_type):
        self.run_type = run_type

    def advance_timestep(
        self,
        nodes,
        time_epoch,
        hour,
        minute,
        second,
        jd,
        jc,
    ):
        """
        Advance the model by one timestep for all stream nodes.

        Apply the applicable StreamNode calculations to
        every node at the provided simulation time.
        """
        t = float(time_epoch)

        if self.run_type == "temperature":
            for node in nodes:
                node.CalcDischarge(t)
            for node in nodes:
                node.CalcHeat(t, hour, minute, second, jd, jc)
            for node in nodes:
                node.maccormick2(t)
            return

        if self.run_type == "solar":
            for node in nodes:
                node.CalcHeat(t, hour, minute, second, jd, jc, True)
            return

        if self.run_type == "hydraulics":
            for node in nodes:
                node.CalcDischarge(t)
            return
