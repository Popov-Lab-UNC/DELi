from dnaio import SequenceRecord


class FailedDecodeAttempt:
    def __init__(self, sequence: SequenceRecord, reason: str):
        self.sequence = sequence
        self.reason = reason

    def __repr__(self):
        return f"{self.__class__.__name__}"
