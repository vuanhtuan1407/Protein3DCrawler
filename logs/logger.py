import logging
import os


class Logger(object):
    def __init__(self, filename=None):
        self.log_dir = os.path.dirname(__file__)
        self.filename = filename
        if filename is not None and not os.path.exists(os.path.join(self.log_dir, filename)):
            logging.warning("Logging directory %s does not exist!", filename)
            self.filename = None
        self.logger = logging.getLogger(os.path.basename(__file__))
        self._configure_logging()

    def _configure_logging(self):
        logging.basicConfig(
            format='%(asctime)s (%(name)s:%(lineno)d) [%(levelname)s]: %(message)s',
            datefmt='%m/%d/%Y %I:%M:%S %p',
            filename=self.filename,
            encoding='utf-8',
            level=logging.DEBUG
        )

    def info(self, msg):
        return self.logger.info(msg)

    def warning(self, msg):
        return self.logger.warning(msg)

    def error(self, msg):
        return self.logger.error(msg)


if __name__ == '__main__':
    logger = Logger()
    logger.info('Hello World!')
