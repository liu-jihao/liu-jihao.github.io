"""Exact documentary decoding and elementary certificate checks."""
from fractions import Fraction
from pathlib import Path
import re
Q = Fraction
ENGINE = 'fractions.Fraction'

def require(condition, message):
    if not condition:
        raise ValueError(message)

def unwrap(text):
    """Physical wraps carry no data; record starts remain on column zero."""
    return re.sub(r'\n[ \t]+', '', text).replace('\r', '')

def record_lines(text):
    return [s.strip() for s in unwrap(text).splitlines() if s.strip()]

def pairs(flat):
    require(len(flat) % 2 == 0, 'odd box length')
    return list(zip(flat[::2], flat[1::2]))

def split_box(box, axis, midpoint):
    require(isinstance(axis, int) and 0 <= axis < len(box)//2, 'invalid split axis')
    lo,hi=box[2*axis:2*axis+2]
    require(lo < midpoint < hi and midpoint == (lo+hi)/2, 'invalid bisection')
    lower,upper=list(box),list(box)
    lower[2*axis+1]=midpoint;upper[2*axis]=midpoint
    return tuple(lower),tuple(upper)

def unique_insert(mapping, key, value):
    require(key not in mapping, 'duplicate record '+str(key))
    mapping[key]=value
