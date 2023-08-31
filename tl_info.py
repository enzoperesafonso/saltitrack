# Define the custom TLInfo class to hold hour angle and track length info of a
# specified declination
class _TLInfo:
    def __init__(self, hour_angle, track_length):
        self.hour_angle = hour_angle
        self.track_length = track_length