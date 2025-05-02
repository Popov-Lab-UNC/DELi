"""base classes for decoding calls"""

import abc


class Call(abc.ABC):
    """base class for all calls"""

    @abc.abstractmethod
    def is_valid(self) -> bool:
        """Return True if call is valid else false"""
        raise NotImplementedError


class ValidCall(Call):
    """base class for all valid calls"""

    def is_valid(self) -> bool:
        """Valid calls are always True"""
        return True


class FailedCall(Call):
    """base class for all failed calls"""

    def is_valid(self) -> bool:
        """Failed calls are always False"""
        return False
