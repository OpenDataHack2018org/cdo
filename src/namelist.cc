/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2018 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "namelist.h"

// Allocates a fresh unused token from the token pull.
static NamelistToken *
namelist_alloc_token(NamelistParser *parser)
{
  const unsigned int TOK_MEM_INCR = 64;

  if (parser->toknext >= parser->num_tokens)
    {
      parser->num_tokens += TOK_MEM_INCR;
      parser->tokens.resize(parser->num_tokens);
    }

  NamelistToken *tok = &parser->tokens[parser->toknext++];
  tok->start = tok->end = -1;
  return tok;
}

// Fills token type and boundaries.
static void
namelist_fill_token(NamelistToken *token, NamelistType type, int start, int end)
{
  token->type = type;
  token->start = start;
  token->end = end;
}

void
namelist_new_object(NamelistParser *parser)
{
  NamelistToken *token;
  token = namelist_alloc_token(parser);
  token->type = NamelistType::OBJECT;
  token->start = parser->pos;
}

// Fills next available token with NAMELIST word.
static NamelistError
namelist_parse_word(NamelistParser *parser, const char *buf, size_t len)
{
  NamelistToken *token;
  int start = parser->pos;

  for (; parser->pos < len && buf[parser->pos] != '\0'; parser->pos++)
    {
      switch (buf[parser->pos])
        {
        case ':':
        case '=':
        case ',':
        case '&':
        case '/':
        case '\r':
        case '\n':
        case '\t':
        case ' ': goto found;
        }

      if (buf[parser->pos] < 32 || buf[parser->pos] >= 127)
        {
          parser->pos = start;
          return NamelistError::INVAL;
        }
    }

found:

  token = namelist_alloc_token(parser);
  namelist_fill_token(token, NamelistType::WORD, start, parser->pos);
  parser->pos--;

  return NamelistError::UNDEFINED;
}

// Fills next token with NAMELIST string.
static NamelistError
namelist_parse_string(NamelistParser *parser, const char *buf, size_t len, char quote)
{
  int start = parser->pos;

  parser->pos++;

  /* Skip starting quote */
  for (; parser->pos < len && buf[parser->pos] != '\0'; parser->pos++)
    {
      char c = buf[parser->pos];

      /* Quote: end of string */
      if (c == quote)
        {
          NamelistToken *token = namelist_alloc_token(parser);
          namelist_fill_token(token, NamelistType::STRING, start + 1, parser->pos);
          return NamelistError::UNDEFINED;
        }

      /* Backslash: Quoted symbol expected */
      if (c == '\\' && parser->pos + 1 < len)
        {
          parser->pos++;
          switch (buf[parser->pos])
            {
            // Allowed escaped symbols
            case '\"':
            case '\\':
            case 'b':
            case 'f':
            case 'r':
            case 'n':
            case 't':
              break;
            // Allows escaped symbol \uXXXX
            case 'u':
              parser->pos++;
              for (int i = 0; i < 4 && parser->pos < len && buf[parser->pos] != '\0'; i++)
                {
                  // If it isn't a hex character we have an error
                  if (!((buf[parser->pos] >= 48 && buf[parser->pos] <= 57) ||  // 0-9
                        (buf[parser->pos] >= 65 && buf[parser->pos] <= 70) ||  // A-F
                        (buf[parser->pos] >= 97 && buf[parser->pos] <= 102)))  // a-f
                    {
                      return NamelistError::INVAL;
                    }
                  parser->pos++;
                }
              parser->pos--;
              break;
            // Unexpected symbol
            default: return NamelistError::INVAL;
            }
        }
    }

  parser->pos = start;

  return NamelistError::PART;
}

static NamelistError
namelist_check_keyname(const char *buf, NamelistToken *t)
{
  switch (t->type)
    {
    case NamelistType::STRING:
      while (isspace((int) buf[t->start]) && t->start < t->end) t->start++;
      while (isspace((int) buf[t->end - 1]) && t->start < t->end) t->end--;
      if ((t->end - t->start) < 1) return NamelistError::EMKEY;
      for (int i = t->start; i < t->end; ++i)
        if (isspace((int) buf[i])) return NamelistError::INKEY;
    case NamelistType::WORD: t->type = NamelistType::KEY; break;
    default: return NamelistError::INTYP;
    }

  return NamelistError::UNDEFINED;
}

NamelistError
NamelistParser::parse(const char *buf, size_t len)
{
  NamelistError status = NamelistError::UNDEFINED;
  NamelistToken *token;

  this->lineno = 1;

  for (; this->pos < len && buf[this->pos] != '\0'; this->pos++)
    {
      char c = buf[this->pos];
      switch (c)
        {
        case '&': namelist_new_object(this); break;
        case '/':
          for (int i = this->toknext - 1; i >= 0; i--)
            {
              token = &this->tokens[i];
              if (token->start != -1 && token->end == -1)
                {
                  if (token->type != NamelistType::OBJECT) return NamelistError::INOBJ;
                  token->end = this->pos + 1;
                  break;
                }
            }
          break;
        case '\t':
        case ' ': break;
        case '\r':
          if (this->pos + 1 < len && buf[this->pos + 1] == '\n') this->pos++;
        case '\n': this->lineno++; break;
        case ',': break;
        case '#':
        case '!':  // Skip to end of line
          for (; this->pos < len && buf[this->pos] != '\0'; this->pos++)
            if (buf[this->pos] == '\r' || buf[this->pos] == '\n')
              {
                this->pos--;
                break;
              }
          break;
        case ':':
        case '=': status = namelist_check_keyname(buf, &this->tokens[this->toknext - 1]); break;
        case '\"':
        case '\'': status = namelist_parse_string(this, buf, len, c); break;
        default: status = namelist_parse_word(this, buf, len); break;
        }

      if (status!=NamelistError::UNDEFINED) return status;
    }

  return status;
}

void
NamelistParser::dump(const char *buf)
{
  unsigned int ntok = this->toknext;
  printf("Number of tokens %d\n", ntok);

  for (unsigned int it = 0; it < ntok; ++it)
    {
      NamelistToken *t = &this->tokens[it];
      int length = t->end - t->start;
      const char *start = buf + t->start;
      printf("Token %u", it + 1);
      if (t->type == NamelistType::OBJECT)
        {
          printf(" NAMELIST=");
          if (length > 80) length = 80;
          printf("'%.*s'", length, start);
        }
      else if (t->type == NamelistType::KEY)
        {
          printf(" KEY=");
          printf("'%.*s'", length, start);
        }
      else if (t->type == NamelistType::WORD)
        {
          printf(" WORD=");
          printf("'%.*s'", length, start);
        }
      else if (t->type == NamelistType::STRING)
        {
          printf(" STRING=");
          printf("'%.*s'", length, start);
        }
      printf("\n");
    }
}

int
NamelistParser::verify()
{
  unsigned int ntok = this->toknext;

  if (ntok)
    {
      NamelistToken *t = &this->tokens[0];
      if (t->type != NamelistType::OBJECT && t->type != NamelistType::KEY) return -1;
    }

  return 0;
}
