"""
Global quiet mode module to handle console output suppression.

This is imported early in each module to ensure consistent quiet mode behavior
across the entire codebase.
"""

import os
import sys
import logging
from typing import Optional

# Global flag to indicate if we're running in quiet mode
QUIET_MODE = False 

# Store original stdout/stderr
ORIGINAL_STDOUT = sys.stdout
ORIGINAL_STDERR = sys.stderr

# Null device for output redirection
NULL_OUTPUT = open(os.devnull, 'w')

# Prevent all logging modules from setting up console loggers in quiet mode
class NullHandler(logging.Handler):
    """A handler that does nothing"""
    def emit(self, record):
        pass

def setup_quiet_mode():
    """Setup quiet mode by silencing all console output from logging"""
    global QUIET_MODE
    
    # Check command line for --quiet flag
    if '--quiet' in sys.argv:
        QUIET_MODE = True
        
        # Disable logging to console
        logging.basicConfig(level=logging.CRITICAL)  # Set highest level to suppress most logs
        root_logger = logging.getLogger()
        
        # Remove any handlers already attached
        for handler in root_logger.handlers[:]:
            root_logger.removeHandler(handler)
            
        # Add a null handler to prevent console output
        root_logger.addHandler(NullHandler())
        
        # Set the global python logging level high to suppress most logs
        null_handler = NullHandler()
        
        # Install null handler on all existing loggers
        for name in logging.root.manager.loggerDict:
            logger = logging.getLogger(name)
            for hdlr in logger.handlers[:]:
                if isinstance(hdlr, logging.StreamHandler) and not isinstance(hdlr, logging.FileHandler):
                    logger.removeHandler(hdlr)
            logger.addHandler(null_handler)

# Run the setup during module import
setup_quiet_mode()