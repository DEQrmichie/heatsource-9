

import logging

from heatsource9.domain.clock import pretty_time


logger = logging.getLogger(__name__)


class ModelRoutine:
    """
    Class that stores the run type and advances all stream nodes by one timestep 
    using the appropriate model run calculations.
    """

    def __init__(self, run_type):
        self.run_type = run_type

    def _warn_dry_channel(self, nodes, time_epoch):
        dry_nodes = [node for node in nodes if node.Q is not None and node.Q < 0.003]
        if not dry_nodes:
            return

        min_q_node = min(dry_nodes, key=lambda node: node.Q)
        logger.warning(
            "The channel is going dry at %s: %d nodes below 0.003 cms, min Q=%.1e cms at km %.2f.",
            pretty_time(time_epoch),
            len(dry_nodes),
            min_q_node.Q,
            min_q_node.km,
        )

    def advance_timestep(
        self,
        nodes,
        time_epoch,
        hour,
        minute,
        second,
        doy,
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
            self._warn_dry_channel(nodes, t)
            for node in nodes:
                node.CalcHeat(t, hour, minute, second, doy, jc)
            for node in nodes:
                node.maccormick2(t)
            return

        if self.run_type == "solar":
            for node in nodes:
                node.CalcHeat(t, hour, minute, second, doy, jc, True)
            return

        if self.run_type == "hydraulics":
            for node in nodes:
                node.CalcDischarge(t)
            self._warn_dry_channel(nodes, t)
            return
