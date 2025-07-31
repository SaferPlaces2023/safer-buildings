import logging

logging.basicConfig(
    format="[%(asctime)s] [%(levelname)-8s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
Logger = logging.getLogger(__name__)
Logger.setLevel(logging.CRITICAL)


def is_debug_mode():
    """
    is_debug_mode - check if the logger is in debug mode
    """
    return Logger.level == logging.DEBUG