from astropy import coordinates



def get_max_track_time(ob : coordinates.HADec):
    """Determines the maximum track length for a target with the specified declination.
    The track length is given as a Universal Time (rather than sidereal time)
    difference. In other words, it is the time difference you measure with a normal
    clock.

    Parameters
    ----------
    declination : float
        Declination of target (in degrees).


    Raises
    ------
    OutOfVisibilityRangeException
        Raised if the target with the given
        declination never enters the visibility zone.

    Returns
    -------
    max_track_length : float
        The maximum tracklength (in seconds).
    """
    pass


def is_visible_from_salt():
    pass