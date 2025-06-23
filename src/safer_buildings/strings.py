# coding=utf-8
# -------------------------------------------------------------------------------------------------
# Copyright (c) 2018-2023
# Developed by Gecosistema.com and Stefano Bagli,Valerio Luzzi for the Gecosistema s.r.l. Agency
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by the Free Software Foundation,
# either version 2 of the License, or (at you option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PORPOSE. See the GNU Gene-
# ral Public License for more details.
# You should have received a copy of the GNU General Public License along with this program. If not,
# see http://www.gnu.org/licenses/.
# -------------------------------------------------------------------------------------------------
import inspect


def isInteger(s):
    """
    isInteger
    """
    try:
        int(s)
        return True
    except ValueError:
        return False


def isFloat(s):
    """
    isFloat
    """
    try:
        float(s)
        return True
    except ValueError:
        return False


def get_default_values(func):
    """
    get_default_values
    """
    signature = inspect.signature(func)
    return {
        k:  v.default if v.default is not inspect.Parameter.empty else None
        for k, v in signature.parameters.items()
    }
    

def parse_event(event, func):
    """
    parse_event
    """
    # Copy the event to avoid side effects
    #kwargs = defaults.copy()
    kwargs = get_default_values(func)


    # Update kwargs with event
    for key in event:
        if key in kwargs:
            kwargs[key] = event[key]    
        else:
            print(f"Option <{key}> is not available")

    # patch numeric and boolean params because can be passed as string
    for key in kwargs:
        value = kwargs[key]
        if isinstance(value, str):
            if value.lower() == "true":
                kwargs[key] = True
            elif value.lower() == "false":
                kwargs[key] = False
            elif isInteger(value):
                kwargs[key] = int(value)
            elif isFloat(value):
                kwargs[key] = float(value)

    return kwargs