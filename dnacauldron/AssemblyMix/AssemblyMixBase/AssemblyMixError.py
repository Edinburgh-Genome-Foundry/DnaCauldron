class AssemblyMixError(Exception):
    
    def __init__(self, message, mix):
        super().__init__(message)
        self.mix = mix