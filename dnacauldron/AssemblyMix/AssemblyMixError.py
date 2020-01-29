class AssemblyMixError(Exception):
    """Simple class to raise and recognize errors happening in Mix operations.

    These have a mix attribute!
    """
    
    def __init__(self, message, mix):
        super().__init__(message)
        self.mix = mix