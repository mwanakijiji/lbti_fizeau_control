def info(text):
    """
    Prints text with a green bold INFO in front of it. Used to print information.

    Usage:
    >> info('This is some information text.')
    """

    color = '\033[34;1m'
    bold = '\033[39;1m'
    normal = '\033[0m'

    print(color + 'INFO: ' + bold + text + normal)


def request(text):
    """
    Prints text with a yellow bold REQUEST in front of it. Used to print requests.

    Usage:
    >> request('This is some request text.')
    """

    color = '\033[31;1m'
    bold = '\033[39;1m'
    normal = '\033[0m'

    print(color + 'REQUEST: ' + bold + text + normal)


