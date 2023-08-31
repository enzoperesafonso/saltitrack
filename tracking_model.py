""" This module is a python implimentation by Enzo Peres Afonso
of the SALT astro ops TrackingModel.java.

A model for tracking with SALT.

 Due to its fixed elevation, SALT can only observe targets
 in a certain visibility zone. Depending on the target's
 declination, the target msy

 a) Cross the visibility zone when rising, leave the zone
     and cross it again when setting. In this case we call
     the first crossing "east track", and the second "west
     track".

  b) Enter the visibility zone when rising, remain in there
     and leave it when setting. In this case we call the
     crossing of the zone "west track".

  c) Not enter the visibility zone at all.

 The east and west track don't depend on the right ascension,
 although it depends on the right ascension and the time of
 the year at what time a target is on the east or west track.

 The tracks for some declination can be obtained from the
 eastTrack(double) and the westTrack(double) method.
 Both methods return hour angles for the start and
 end of the track.

 When actually observing a target, SALT's tracker must follow
 the target, and this can only be done for a limited time,
 the track length. The track length depends not only on the
 target's declination, but on its (local) hour angle as well.
 In general it is shorter than the length of the east or west
 track.

 The track length can be obtained by means of the
 trackLength(double, double)} method.
 Note that there are two slightly different meanings of the
 word "track". In "east track" or "west track" it denotes the
 time the target spends in the visibility zone. In "track
 length" it denotes the time while the target can be observed. """


import bisect
import exceptions



# Define the custom TLInfo class to hold hour angle and track length info of a
# specified declination
class _TLInfo:
    def __init__(self, hour_angle, track_length):
        self.hour_angle = hour_angle
        self.track_length = track_length


# Define the custom TrackHA class to hold hour angle at start and end of a
# specified track
class _TrackHA:
    def __init__(self, start, end):
        self.start = start
        self.end = end


# Constants
_UT_TO_ST_FACTOR = 1.0027379  # Conversion factor from UT to sidereal time


class TrackingModel:
    """    This class is a python implimentation of the SALT astro ops
    TrackingModel.java. A model for tracking with SALT.

    Attributes
    ----------
    declinations : float
        An ordered list of declinations observable by SALT.
    east_track : _TLInfo
        A list of TLInfo objects 
    west_tracks : _TLInfo
        A list of TLInfo objects 

    Methods     # DONT INCLUDE PRIVATE METHODS IN DOC STRINGs
    -------
    read_model():
        Reads and imports the model data from visDataOrdered.dat.

    lower_index_for_declination(declination: float):
        Determines the maxmimum declination index i with {declination[i] ≤ declination}.

    east_track(declination: float):
        Determines the start and end hour angles of the east track for the specified
        target declination.

    west_track(declination: float):
        Determines the start and end hour angles of the the west track for the specified
        target declination.
    """



    def __init__(self):
        self.declinations = []
        self.east_tracks = []
        self.west_tracks = []
        self.read_model()

    def read_model(self):
        # Read in the columns from the data file.
        with open("visDataOrdered.dat", "r") as file:
            # Skip the header line
            next(file)

            column1 = []  # declination
            column2 = []  # local hour angle
            column3 = []  # track length

            for line in file:
                declination, hour_angle, track_length, azimuth_angle = map(
                    float, line.strip().split()
                )

                column1.append(declination)
                column2.append(hour_angle)
                column3.append(track_length)

            # Process the data and populate east_tracks and west_tracks lists
            east_track_info = None
            west_track_info = None
            previous_declination = None

            for i in range(len(column1)):
                declination = column1[i]
                if previous_declination is None or declination != previous_declination:
                    # next declination.
                    self.declinations.append(declination)
                    east_track_info = []
                    self.east_tracks.append(east_track_info)
                    west_track_info = []
                    self.west_tracks.append(west_track_info)
                previous_declination = declination

                # If the target has a declination between -62.75 deg and -1.75 deg, there are two distinct tracks.
                # Otherwise only an east track.
                hour_angle = column2[i]
                track_length = column3[i]
                tl_info = _TLInfo(hour_angle, track_length)
                if -62.75 <= declination < -1.75 and hour_angle > 0:
                    west_track_info.append(tl_info)
                else:
                    east_track_info.append(tl_info)

            # Let a, b be the hour angle and track length in the last line for a declination in the data file.
            # Then a + b is the end of the track. We explicitly add this to our data.
            for tl_info_list in self.east_tracks:
                last_tl_info = tl_info_list[-1]
                end_of_track = (
                    last_tl_info.hour_angle
                    + _UT_TO_ST_FACTOR * last_tl_info.track_length / 3600.0
                )
                tl_info_list.append(_TLInfo(end_of_track, 0.0))

            for tl_info_list in self.west_tracks:
                # Only add the track end if there actually is a west track.
                if tl_info_list:
                    last_tl_info = tl_info_list[-1]
                    end_of_track = (
                        last_tl_info.hour_angle
                        + _UT_TO_ST_FACTOR * last_tl_info.track_length / 3600.0
                    )
                    tl_info_list.append(_TLInfo(end_of_track, 0.0))

    def lower_index_for_declination(self, declination: float):     ##### PLACE TO CONVERT ANGLE TO FLOAT
        """Determines the maxmimum declination index i with {declinations[i] ≤ declination}.

        Parameters
        ----------
        value : float
            value (float): Value to get closest index of


        Raises
        ------
        OutOfVisibilityRangeException
            Raised if the declination is out of the bounds of the models declinations.

        Returns
        -------
        index : int
            The lowest index of the value in the ordered list of declinations
        """

        if declination < self.declinations[0] or declination > self.declinations[-1]:
            message = f"The value {declination} lies outside the allowed range from {self.declinations[0]} to {self.declinations[-1]}"
            raise exceptions.OutOfVisibilityRangeException(message)

        i = bisect.bisect_right(self.declinations, declination) - 1
        if i == len(self.declinations) - 1:
            i -= 1
        return i

    def east_track(self, declination: float):
        """Determines the start and end hour angles of the east track for the specified
        target declination.Throws IllegalArgumentException if the target with the given
        declination never enters the visibility zone.

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
        TrackHA : _TrackHA (or subclass)
            Object storing the start and end hour angles of the east track for the
            specified declination.
        """
        lower_index = self.lower_index_for_declination(declination)   
        east_track1 = self.east_tracks[lower_index]
        east_track2 = self.east_tracks[lower_index + 1]
        west_track1 = self.west_tracks[lower_index]
        west_track2 = self.west_tracks[lower_index + 1]
        declination1 = self.declinations[lower_index]
        declination2 = self.declinations[lower_index + 1]

        if west_track1 and west_track2:
            # There are west tracks throughout the declination interval.
            start1 = east_track1[0].hour_angle
            end1 = east_track1[-1].hour_angle
            start2 = east_track2[0].hour_angle
            end2 = east_track2[-1].hour_angle
        elif not west_track1 and west_track2:
            # We assume that at the lower index the east and west
            # track touch each other (at the hour angle 0).
            start1 = east_track1[0].hour_angle
            end1 = 0
            start2 = east_track2[0].hour_angle
            end2 = east_track2[-1].hour_angle
        elif west_track1 and not west_track2:
            # We assume that at the upper index the east and west
            # track touch each other (at the hour angle 0).
            start1 = east_track1[0].hour_angle
            end1 = east_track1[-1].hour_angle
            start2 = east_track2[0].hour_angle
            end2 = 0
        else:
            # There are no west tracks throughout the declination interval.
            start1 = east_track1[0].hour_angle
            end1 = east_track1[-1].hour_angle
            start2 = east_track2[0].hour_angle
            end2 = east_track2[-1].hour_angle

        # Perform a linear interpolation.
        start = start1 + ((start2 - start1) / (declination2 - declination1)) * (
            declination - declination1
        )
        end = end1 + ((end2 - end1) / (declination2 - declination1)) * (
            declination - declination1
        )

        return _TrackHA(start, end)

    def west_track(self, declination: float):
        """Determines the start and end hour angles of the the west track for the specified
        target declination. If there is no west track, None is returned.

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
        TrackHA : _TrackHA (or subclass)
            Object storing the start and end hour angles of the east track for the
            specified declination.
        """

        lower_index = self.lower_index_for_declination(declination)
        east_track1 = self.east_tracks[lower_index]
        east_track2 = self.east_tracks[lower_index + 1]
        west_track1 = self.west_tracks[lower_index]
        west_track2 = self.west_tracks[lower_index + 1]
        declination1 = self.declinations[lower_index]
        declination2 = self.declinations[lower_index + 1]

        if west_track1 and west_track2:
            # There are west tracks throughout the declination interval.
            start1 = west_track1[0].hour_angle
            end1 = west_track1[-1].hour_angle
            start2 = west_track2[0].hour_angle
            end2 = west_track2[-1].hour_angle
        elif not west_track1 and west_track2:
            # We assume that at the lower index the east and west
            # track touch each other (at the hour angle 0).
            start1 = 0
            end1 = east_track1[-1].hour_angle
            start2 = west_track2[0].hour_angle
            end2 = west_track2[-1].hour_angle
        elif west_track1 and not west_track2:
            # We assume that at the upper index the east and west
            # track touch each other (at the hour angle 0).
            start1 = west_track1[0].hour_angle
            end1 = west_track1[-1].hour_angle
            start2 = 0
            end2 = east_track2[-1].hour_angle
        else:
            # There are no west tracks throughout the declination interval.
            return None

        # Perform a linear interpolation.
        start = start1 + ((start2 - start1) / (declination2 - declination1)) * (
            declination - declination1
        )
        end = end1 + ((end2 - end1) / (declination2 - declination1)) * (
            declination - declination1
        )

        return _TrackHA(start, end)

    def track_length_for_hour_angle(self, hour_angle: float, tlInfoList: _TLInfo):
        """Determines the track length from a _TLInfo object for a given hour angle.

        Parameters
        ----------
        hour_angle : float
            The hour angle for which the track length is to be determined.
        tlInfoList : _TLInfo
            The track length object of a specific declination for which the track
            length for hour_angle will be determined.

        Returns
        -------
        track_length : float
            The track length for the given hour_angle from the given _TLInfo objcet.
        """

        if (
            hour_angle < tlInfoList[0].hour_angle
            or hour_angle >= tlInfoList[-1].hour_angle
        ):
            return 0

        # Find the index hour_angle_index so that HA(hour_angle_index)
        # <= hour_angle <= HA(hour_angle_index + 1).
        hour_angle_index = 0
        while tlInfoList[hour_angle_index].hour_angle <= hour_angle:
            hour_angle_index += 1
        hour_angle_index -= 1
        if hour_angle_index == len(tlInfoList) - 1:
            hour_angle_index -= 1

        # Find the track length by linear interpolation.
        hour_angle_1 = tlInfoList[hour_angle_index].hour_angle
        track_length_1 = tlInfoList[hour_angle_index].track_length
        hour_angle_2 = tlInfoList[hour_angle_index + 1].hour_angle
        track_length_2 = tlInfoList[hour_angle_index + 1].track_length

        return track_length_1 + (
            (track_length_2 - track_length_1) / (hour_angle_2 - hour_angle_1)
        ) * (hour_angle - hour_angle_1)

    def track_length(self, declination: float, hour_angle: float) -> float:
        """Determines the track length for a target with the specified declination and
        specified hour angle. The track length is given as a Universal Time
        (rather than sidereal time) difference. In other words, it is the time
        difference you measure with a normal clock.

        Parameters
        ----------
        declination : float
            Declination of target (in degrees).

        hour_angle : float
            Hour angle of target (as a fractional hour value).


        Raises
        ------
        OutOfVisibilityRangeException
            Raised if the target with the given
            declination never enters the visibility zone.

        Returns
        -------
        track_length : float
            The tracklength (in seconds).
        """

        declinationIndex = self.lower_index_for_declination(declination)
        
        eastTrackInfoList1 = self.east_tracks[declinationIndex]
        westTrackInfoList1 = self.west_tracks[declinationIndex]
        eastTrackInfoList2 = self.east_tracks[declinationIndex + 1]
        westTrackInfoList2 = self.west_tracks[declinationIndex + 1]

        eastTrack = self.east_track(declination)
        westTrack = self.west_track(declination)

        if westTrack is not None:
            if (
                hour_angle < eastTrack.start
                or (hour_angle > eastTrack.end and hour_angle < westTrack.start)
                or hour_angle > westTrack.end
            ):
                return 0
        else:
            if hour_angle < eastTrack.start or hour_angle > eastTrack.end:
                return 0

        if westTrackInfoList1 and westTrackInfoList2:
            tlInfoList1 = eastTrackInfoList1 if hour_angle < 0 else westTrackInfoList1
            tlInfoList2 = eastTrackInfoList2 if hour_angle < 0 else westTrackInfoList2
            track_length_1 = self.track_length_for_hour_angle(hour_angle, tlInfoList1)
            track_length_2 = self.track_length_for_hour_angle(hour_angle, tlInfoList2)
        elif not westTrackInfoList1 and westTrackInfoList2:
            track_length_1 = self.track_length_for_hour_angle(
                hour_angle, eastTrackInfoList1
            )
            if hour_angle < 0:
                remaining_east_track_time = 3600 * (0 - hour_angle) / _UT_TO_ST_FACTOR
                if remaining_east_track_time < track_length_1:
                    track_length_1 = remaining_east_track_time
            tlInfoList2 = eastTrackInfoList2 if hour_angle < 0 else westTrackInfoList2
            track_length_2 = self.track_length_for_hour_angle(hour_angle, tlInfoList2)
        elif westTrackInfoList1 and not westTrackInfoList2:
            tlInfoList1 = eastTrackInfoList1 if hour_angle < 0 else westTrackInfoList1
            track_length_1 = self.track_length_for_hour_angle(hour_angle, tlInfoList1)
            track_length_2 = self.track_length_for_hour_angle(
                hour_angle, eastTrackInfoList2
            )
            if hour_angle < 0:
                remaining_east_track_time = 3600 * (0 - hour_angle) / _UT_TO_ST_FACTOR
                if remaining_east_track_time < track_length_2:
                    track_length_2 = remaining_east_track_time
        else:
            track_length_1 = self.track_length_for_hour_angle(
                hour_angle, eastTrackInfoList1
            )
            track_length_2 = self.track_length_for_hour_angle(
                hour_angle, eastTrackInfoList2
            )

        declination_1 = self.declinations[declinationIndex]
        declination_2 = self.declinations[declinationIndex + 1]
        return track_length_1 + (
            (track_length_2 - track_length_1) / (declination_2 - declination_1)
        ) * (declination - declination_1)

    def maximum_track_length(self, declination: float):
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


        declination_index = self.lower_index_for_declination( 
                declination)
            

        # Collect the relevant track information
        declination_index = self.lower_index_for_declination(declination)
        tlInfoList = []
        tlInfoList.extend(self.east_tracks[declination_index])
        tlInfoList.extend(self.west_tracks[declination_index])
        tlInfoList.extend(self.east_tracks[declination_index + 1])
        tlInfoList.extend(self.west_tracks[declination_index + 1])

        # Extract the hour angles
        hour_angles = [tlInfo.hour_angle for tlInfo in tlInfoList]

        # Calculate the maximum track length for these hour angles
        maximum_track_length = 0
        for hour_angle in hour_angles:
            track_length = self.track_length(declination, hour_angle)
            if track_length > maximum_track_length:
                maximum_track_length = track_length

        # As we are using linear interpolation in the track
        # length calculation, this maximum should be as good
        # a value as we can get
        return maximum_track_length
