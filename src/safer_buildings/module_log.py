import logging

logging.basicConfig(format="[%(levelname)-8s] %(message)s")
Logger = logging.getLogger(__name__)
Logger.setLevel(logging.CRITICAL)


def is_debug_mode():
    """
    is_debug_mode - check if the logger is in debug mode
    """
    return Logger.level == logging.DEBUG