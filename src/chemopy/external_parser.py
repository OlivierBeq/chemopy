# -*- coding: utf-8 -*-

"""Definition of external.cfg file parser for external tools."""

import os
from configparser import ConfigParser, SectionProxy
from typing import List


class ExternalToolsParser:
    """A parser for the config file of external tools."""

    def __init__(self, path: str = None,
                 required_fields: List[str] = None,
                 skip_errors: bool = False) -> None:
        """Initialize an ExternalToolsParser.

        :param path: path to configuration file. If None,
                     the default config file, distributed
                     with chemopy is used in place.
        """
        if path is None:
            path = os.path.join(os.path.dirname(__file__), 'external.cfg')
        cfp = ConfigParser()
        if not cfp.read(path):
            raise FileNotFoundError(f'Config file not found: {path}')
        test = ExternalToolsParser._test_required_fieds_of_section
        for section in cfp.sections():
            if required_fields and not (test(cfp[section], required_fields) or skip_errors):
                raise AttributeError(f'Required attributes {required_fields} are not defined in {path}.')
        self.tools = cfp._sections
        self.path = path

    def update(self, path: str = None) -> None:
        """Update the parser with new sections and values since last read.

        :param path: path to configuration file. If None,
                     the previous config file, is reread.
        """
        self.__init__(self.path)

    @staticmethod
    def _test_required_fieds_of_section(section: SectionProxy, required_fields: List[str]) -> bool:
        """Test that all required fields are present in all sections of the parsed ConfigParser.

        :param section: section of a ConfigPArser to test the fields of
        :param required_fields: fields required in all sections of the ConfigParser
        """
        return len(set(section.keys()) & set(required_fields)) >= len(required_fields)
