# Licensed under a 3-clause BSD style license - see LICENSE


import html
import re

# import astropy.io.ascii as asciitable
from astropy.table import Table

from astroquery.query import BaseQuery
from astroquery.utils import async_to_sync, prepend_docstr_nosections
from . import conf
from astroquery.exceptions import TableParseError

__all__ = ['NistLevels', 'NistLevelsClass']


def _strip_blanks(table):
    """
    Remove blank lines from table (included for "human readability" but
    useless to us...
    returns a single string joined by \n newlines

    Parameters
    ----------
    table : str
       table to strip as a string

    Returns
    -------
    single string joined by newlines.
    """
    numbersletters = re.compile("[0-9A-Za-z]")
    if isinstance(table, str):
        table = table.split('\n')
    table = [line for line in table if numbersletters.search(line)]
    return "\n".join(table)


@async_to_sync
class NistLevelsClass(BaseQuery):
    URL = conf.server
    TIMEOUT = conf.timeout
    unit_code = {'Angstrom': 0,
                 'nm': 1,
                 'um': 2}
    energy_level_code = {'cm-1': 0, 'invcm': 0, 'cm': 0,
                         'ev': 1, 'eV': 1, 'EV': 1, 'electronvolt': 1,
                         'R': 2, 'Rydberg': 2, 'rydberg': 2}

    def _args_to_payload(self, *args, **kwargs):
        """
        Serves the same purpose as `~NistClass.query` but returns
        the raw HTTP response rather than a `~astropy.table.Table` object.

        Parameters
        ----------
        linename : str, optional
            The spectrum to fetch. Defaults to "H I"
        energy_level_unit : str, optional
            The energy level units must be one of the following:
            'R', 'Rydberg', 'rydberg', 'cm', 'cm-1', 'EV', 'eV',
            'electronvolt', 'ev', 'invcm' Defaults to 'eV'.
        get_query_payload : bool, optional
            If true then returns the dictionary of query parameters, posted to
            remote server. Defaults to `False`.

        Returns
        -------
        request_payload : dict
            The dictionary of parameters sent with the HTTP request

        """
        request_payload = {}
        linename = kwargs["linename"]
        request_payload["de"] = 0
        request_payload["spectrum"] = linename
        request_payload["submit"] = "Retrieve Data"
        request_payload["units"] = NistLevels.energy_level_code[
            kwargs["energy_level_unit"]]
        request_payload["format"] = 1  # ascii
        request_payload["output"] = 0  # entirely rather than pagewise
        request_payload["page_size"] = 15
        request_payload["multiplet_ordered"] = 0
        request_payload["conf_out"] = "on"
        request_payload["term_out"] = "on"
        request_payload["level_out"] = "on"
        request_payload["unc_out"] = "on"
        request_payload["j_out"] = "on"
        request_payload["g_out"] = "on"
        request_payload["lande_out"] = "on"
        request_payload["perc_out"] = "on"
        request_payload["biblio"] = "on"
        request_payload["splitting"] = 1
        return request_payload

    @prepend_docstr_nosections("\n" + _args_to_payload.__doc__)
    def query_async(self, linename="H I", energy_level_unit='eV', get_query_payload=False):
        """
        Returns
        -------
        response : `requests.Response` object
            The response of the HTTP request.
        """
        request_payload = self._args_to_payload(
            linename=linename,
            energy_level_unit=energy_level_unit)
        if get_query_payload:
            return request_payload

        response = self._request("GET", url=NistLevels.URL, params=request_payload,
                                 timeout=NistLevels.TIMEOUT)
        return response

    def _parse_result(self, response, *, verbose=False):
        """
        Parses the results from the HTTP response to `astropy.table.Table`.

        Parameters
        ----------
        response : `requests.Response`
            The HTTP response object

        Returns
        -------
        table : `astropy.table.Table`
        """

        pre_re = re.compile("<PRE>(.*)</PRE>", flags=re.DOTALL)
        links_re = re.compile(r"<\/?a(.|\n)*?>")
        content = str(response.text)

        try:
            pre = pre_re.findall(content)[0]
        except IndexError:
            raise Exception("Result did not contain a table")
        try:
            table = _strip_blanks(pre)
            table = links_re.sub(r'\1', table)
            table = html.unescape(table)
            table = Table.read(table, format='ascii.fixed_width', data_start=1, delimiter='|')

            return table
        except Exception as ex:
            self.response = response
            self.table_parse_error = ex
            raise TableParseError("Failed to parse asciitable! The raw "
                                  "response can be found in self.response, "
                                  "and the error in self.table_parse_error.")


NistLevels = NistLevelsClass()
