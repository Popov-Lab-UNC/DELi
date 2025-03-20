"""base classes for decoding calls"""

import abc


class Call(abc.ABC):
    """base class for all calls"""

    @abc.abstractmethod
    def valid_call(self) -> bool:
        """Return True if call is valid else false"""
        raise NotImplementedError


class ValidCall(Call):
    """base class for all valid calls"""

    def valid_call(self) -> bool:
        """Valid calls are always True"""
        return True


class FailedCall(Call):
    """base class for all failed calls"""

    def valid_call(self) -> bool:
        """Failed calls are always False"""
        return False
