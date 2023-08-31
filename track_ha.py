# Define the custom TrackHA class to hold hour angle at start and end of a
# specified track
class _TrackHA:
    def __init__(self, start, end):
        self.start = start
        self.end = end