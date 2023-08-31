# Raised when the input declination is not reachable by SALT"
class OutOfVisibilityRangeException(Exception):
    def __init__(self, message):
        super().__init__(message)