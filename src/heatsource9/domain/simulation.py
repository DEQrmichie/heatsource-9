from heatsource9.domain.clock import Clock


class Simulation:
    """
    Model simulation container. It stores the clock and node list.

    Simulation is produced by ModelSetup.build() and is used by 
    the run.model_runner. Simulation itself does not perform computation.
    """
    def __init__(self, clock, nodes):
        self.clock = clock
        self.nodes = nodes
