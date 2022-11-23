"""
relaxedjson (https://github.com/simon-engledew/relaxedjson) is a library
for parsing JSON that is missing quotes around its keys.
Uses the parser-combinator library Parsec (https://github.com/sighingnow/parsec.py)
To install from pypi::
    pip install relaxedjson
To install as an egg-link in development mode::
    python setup.py develop -N
"""

import re
from parsec import (
    sepBy,
    regex,
    string,
    generate,
    many,
    endBy
)

whitespace = regex(r'[^\S\r\n]*', re.MULTILINE)
vertwhitespace = regex(r'\s*', re.MULTILINE)

lexeme = lambda p: p << whitespace
openlexeme = lambda p: p << vertwhitespace
closelexeme = lambda p: vertwhitespace >> p << whitespace

comment = string('/*') >> regex(r'(?:[^*]|\*(?!\/))+', re.MULTILINE) << string('*/')
comment = openlexeme(comment)

lbrace = openlexeme(string('{'))
rbrace = closelexeme(string('}'))
lbrack = openlexeme(string('['))
rbrack = closelexeme(string(']'))
lbround = openlexeme(string('('))
rbround = closelexeme(string(')'))
colon = lexeme(string(':'))
separator = openlexeme(string(',')) | vertwhitespace
true = lexeme(string('True')).result(True)
false = lexeme(string('False')).result(False)
null = lexeme(string('None')).result(None)
quote = string('"') | string("'")
literal = lexeme(regex(r'(?!True)(?!False)(?!None)[a-zA-Z][-_a-zA-Z0-9.]*'))

def number_float():
    return lexeme(
        regex(r'-?(0|[1-9][0-9]*)([.][0-9]+)([eE][+-]?[0-9]+)?')
    ).parsecmap(float)

def number_int():
    return lexeme(
        regex(r'-?(0|[1-9][0-9]*)')
    ).parsecmap(int)

def to_chr(value):
    return chr(int(value[1:], 16))

def charseq(end_quote):
    def string_part():
        return regex(r'[^{}\\]+'.format(end_quote))

    def string_esc():
        return string('\\') >> (
            string('\\')
            | string('/')
            | string('b').result('\b')
            | string('f').result('\f')
            | string('n').result('\n')
            | string('r').result('\r')
            | string('t').result('\t')
            | regex(r'u[0-9a-fA-F]{4}').parsecmap(to_chr)
            | string(end_quote)
        )
    return string_part() | string_esc()

class StopGenerator(StopIteration):
    def __init__(self, value):
        self.value = value

@lexeme
@generate
def quoted():
    end_quote = yield quote
    body = yield many(charseq(end_quote))
    yield string(end_quote)
    raise StopGenerator(''.join(body))

@generate
def array():
    yield lbrack << many(comment)
    elements = yield sepBy(value, separator)
    yield rbrack << many(comment)
    raise StopGenerator(elements)

@generate
def tuple_object():
    yield lbround << many(comment)
    elements = yield sepBy(value, separator)
    yield rbround << many(comment)
    raise StopGenerator(tuple(elements))


@generate
def object_pair():
    yield many(comment)
    key = yield quoted | literal
    yield many(comment) << colon << many(comment)
    val = yield value
    yield many(comment)
    raise StopGenerator((key, val))


@generate
def json_object():
    yield lbrace << many(comment)
    pairs = yield sepBy(object_pair, separator)
    yield many(comment) << rbrace
    raise StopGenerator(dict(pairs))

value = literal | quoted | number_float() | number_int() | json_object | array | tuple_object | true | false | null
value = many(comment) >> value << many(comment)

parser = vertwhitespace >> (json_object | array | tuple_object) << vertwhitespace

def parse(text):
    """
    Attempt to parse JSON text, returning an object
    :param text: String containing json
    :rtype: dict or list
    :raises: parsec.ParseError
    """
    return parser.parse_strict(text)

__all__ = ['parse']
