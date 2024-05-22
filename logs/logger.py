import logging
import os


class Logger(object):
    def __init__(self, filename=None, use_console: bool = False):
        # reset logging to default
        logging.root.handlers = []
        self.log_dir = os.path.dirname(__file__)
        self.filename = filename
        self.use_console = use_console
        self.logger = logging.getLogger(os.path.basename(__name__))
        self._configure_logging()

    def _configure_logging(self):
        formatter = logging.Formatter(
            fmt='%(asctime)s (%(name)s:%(lineno)d) [%(levelname)s]: %(message)s',
            datefmt='%m/%d/%Y %I:%M:%S %p'
        )

        file_handler = logging.FileHandler(self.filename, encoding='utf-8', delay=False)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

        if self.use_console:
            stream_handler = logging.StreamHandler()
            stream_handler.setFormatter(formatter)
            self.logger.addHandler(stream_handler)

        self.logger.setLevel(logging.DEBUG)

    def info(self, msg):
        return self.logger.info(msg)

    def warning(self, msg):
        return self.logger.warning(msg)

    def error(self, msg):
        return self.logger.error(msg)


if __name__ == '__main__':
    logger = Logger(filename='./test.log', use_console=True)
    logger.info('Hello World!')
