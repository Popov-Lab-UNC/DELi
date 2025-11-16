import abc

from dnaio import SequenceRecord


class FailedDecodeAttempt:
    def __init__(self, sequence: SequenceRecord, reason: str):
        self.sequence = sequence
        self.reason = reason

    def __repr__(self):
        return f"{self.__class__.__name__}"


class Call(abc.ABC):
    """base class for all calls"""

    @abc.abstractmethod
    def is_valid(self) -> bool:
        """Return True if call is valid else false"""
        raise NotImplementedError
