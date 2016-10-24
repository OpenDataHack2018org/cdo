#include "jsmn.h"

/**
 * Allocates a fresh unused token from the token pull.
 */
static
jsmntok_t *jsmn_alloc_token(jsmn_parser *parser, jsmntok_t *tokens, size_t num_tokens)
{
  if (parser->toknext >= num_tokens) return NULL;

  jsmntok_t *tok = &tokens[parser->toknext++];
  tok->start = tok->end = -1;
  tok->size = 0;
  return tok;
}

const char *Types[] = {"Undefined", "Object", "Array", "String", "Primitive"};

const char *JS;
/**
 * Fills token type and boundaries.
 */
static
void jsmn_fill_token(jsmntok_t *token, jsmntype_t type, int start, int end)
{
  static int tok = 0;
  token->type = type;
  token->start = start;
  token->end = end;
  token->size = 0;
  printf("    token %d: start %d end %d %s >%.*s<\n", ++tok, start, end, Types[type],  end - start, JS+start);
}

/**
 * Fills next available token with JSON primitive.
 */
static int jsmn_parse_primitive(jsmn_parser *parser, const char *js, size_t len, jsmntok_t *tokens, size_t num_tokens)
{
  jsmntok_t *token;
  int start = parser->pos;

  for (; parser->pos < len && js[parser->pos] != '\0'; parser->pos++)
    {
      switch (js[parser->pos])
        {
#ifndef JSMN_STRICT
        /* In strict mode primitive must be followed by "," or "}" or "]" */
        case ':': case '=':
#endif
        case '\t' : case '\r' : case '\n' : case ' ' :
        case ','  : case ']'  : case '}' : case '/' : 
          goto found;
        }
      if (js[parser->pos] < 32 || js[parser->pos] >= 127)
        {
          parser->pos = start;
          return JSMN_ERROR_INVAL;
        }
    }
#ifdef JSMN_STRICT
  /* In strict mode primitive must be followed by a comma/object/array */
  parser->pos = start;
  return JSMN_ERROR_PART;
#endif

 found:
  if (tokens == NULL)
    {
      parser->pos--;
      return 0;
    }
  token = jsmn_alloc_token(parser, tokens, num_tokens);
  if (token == NULL)
    {
      parser->pos = start;
      return JSMN_ERROR_NOMEM;
    }
  jsmn_fill_token(token, JSMN_PRIMITIVE, start, parser->pos);
  parser->pos--;
  return 0;
}

/**
 * Fills next token with JSON string.
 */
static int jsmn_parse_string(jsmn_parser *parser, const char *js, size_t len, jsmntok_t *tokens, size_t num_tokens)
{
  jsmntok_t *token;
  int start = parser->pos;

  parser->pos++;

  /* Skip starting quote */
  for (; parser->pos < len && js[parser->pos] != '\0'; parser->pos++)
    {
      char c = js[parser->pos];

      /* Quote: end of string */
      if (c == '\"')
        {
          if (tokens == NULL) return 0;

          token = jsmn_alloc_token(parser, tokens, num_tokens);
          if (token == NULL)
            {
              parser->pos = start;
              return JSMN_ERROR_NOMEM;
            }
          jsmn_fill_token(token, JSMN_STRING, start+1, parser->pos);
          return 0;
        }

      /* Backslash: Quoted symbol expected */
      if (c == '\\' && parser->pos + 1 < len)
        {
          parser->pos++;
          switch (js[parser->pos]) {
            /* Allowed escaped symbols */
          case '\"': /* case '/' : */ case '\\' : case 'b' :
          case 'f' : case 'r' : case 'n'  : case 't' :
            break;
            /* Allows escaped symbol \uXXXX */
          case 'u':
            parser->pos++;
            for(int i = 0; i < 4 && parser->pos < len && js[parser->pos] != '\0'; i++)
              {
                /* If it isn't a hex character we have an error */
                if(!((js[parser->pos] >= 48 && js[parser->pos] <= 57) || /* 0-9 */
                     (js[parser->pos] >= 65 && js[parser->pos] <= 70) || /* A-F */
                     (js[parser->pos] >= 97 && js[parser->pos] <= 102))) { /* a-f */
                  parser->pos = start;
                  return JSMN_ERROR_INVAL;
                }
                parser->pos++;
              }
            parser->pos--;
            break;
            /* Unexpected symbol */
          default:
            parser->pos = start;
            return JSMN_ERROR_INVAL;
          }
        }
    }
  parser->pos = start;
  return JSMN_ERROR_PART;
}

/**
 * Parse JSON string and fill tokens.
 */
int jsmn_parse(jsmn_parser *parser, const char *js, size_t len, jsmntok_t *tokens, unsigned int num_tokens)
{
  JS=js;
  int r;
  int i;
  jsmntok_t *token;
  int count = parser->toknext;

  if ( parser->pos == 0 )
    {
          count++;
          if (tokens == NULL) return JSMN_ERROR_NOMEM;
          token = jsmn_alloc_token(parser, tokens, num_tokens);
          if (token == NULL) return JSMN_ERROR_NOMEM;
          if (parser->toksuper != -1)
            {
              tokens[parser->toksuper].size++;
            }
          token->type = JSMN_OBJECT;
          token->start = parser->pos;
          printf("xstart object %d\n", token->start);
          parser->toksuper = parser->toknext - 1;
    }

  for (; parser->pos < len && js[parser->pos] != '\0'; parser->pos++)
    {
      jsmntype_t type;

      char c = js[parser->pos];
      switch (c)
        {
        case '&':
          /*
          parser->pos++;
          r = jsmn_parse_primitive(parser, js, len, tokens, num_tokens);
          if (r < 0) return r;
          count++;
          if (parser->toksuper != -1 && tokens != NULL)
            tokens[parser->toksuper].size++;
          */

          // parser->pos++;
          if (tokens == NULL) break;
          // count++;
          //  if (parser->toksuper != -1 && tokens != NULL)
          //  tokens[parser->toksuper].size++;
          token = jsmn_alloc_token(parser, tokens, num_tokens);
          if (token == NULL)
            {
              printf("return1\n");
              //parser->pos = start;
              return JSMN_ERROR_NOMEM;
            }
          jsmn_fill_token(token, JSMN_PRIMITIVE, parser->pos+1, parser->pos+1);
          //parser->pos++;
          count++;
          if (parser->toksuper != -1 && tokens != NULL)
            tokens[parser->toksuper].size++;
 
          count++;
          if (tokens == NULL) break;
          token = jsmn_alloc_token(parser, tokens, num_tokens);
          if (token == NULL)
            {
              return JSMN_ERROR_NOMEM;
            }
          if (parser->toksuper != -1)
            {
              tokens[parser->toksuper].size++;
            }
          token->type = JSMN_OBJECT;
          token->start = parser->pos;
          printf("  start object %d\n", token->start);
          parser->toksuper = parser->toknext - 1;

          break;
        case '/':
          if (tokens == NULL) break;
          type = JSMN_OBJECT;
          for (i = parser->toknext - 1; i >= 0; i--)
            {
              token = &tokens[i];
              if (token->start != -1 && token->end == -1)
                {
                  if (token->type != type) return JSMN_ERROR_INVAL;
                  parser->toksuper = -1;
                  token->end = parser->pos + 1;
                  break;
                }
            }
          printf("  end %d\n", token->end);
          /* Error if unmatched closing bracket */
          if (i == -1) return JSMN_ERROR_INVAL;
          for (; i >= 0; i--)
            {
              token = &tokens[i];
              if (token->start != -1 && token->end == -1) {
                parser->toksuper = i;
                break;
              }
            }
          break;
        case '\"':
          r = jsmn_parse_string(parser, js, len, tokens, num_tokens);
          if (r < 0) return r;
          count++;
          if (parser->toksuper != -1 && tokens != NULL)
            tokens[parser->toksuper].size++;
          break;
        case '#': case '!': // Skip to end of line
          for (; parser->pos < len && js[parser->pos] != '\0'; parser->pos++)
            if ( js[parser->pos] == '\r' || js[parser->pos] == '\n' ) break;
        case '\t' : case '\r' : case '\n' : case ' ':
          break;
        case ':': case '=':
          parser->toksuper = parser->toknext - 1;
          break;
        case ',':
          if (tokens != NULL && parser->toksuper != -1 &&
              tokens[parser->toksuper].type != JSMN_ARRAY &&
              tokens[parser->toksuper].type != JSMN_OBJECT)
            {
              for (i = parser->toknext - 1; i >= 0; i--)
                {
                  if (tokens[i].type == JSMN_ARRAY || tokens[i].type == JSMN_OBJECT)
                    {
                      if (tokens[i].start != -1 && tokens[i].end == -1)
                        {
                          parser->toksuper = i;
                          break;
                        }
                    }
                }
            }
          break;
#ifdef JSMN_STRICT
          /* In strict mode primitives are: numbers and booleans */
        case '-': case '0': case '1' : case '2': case '3' : case '4':
        case '5': case '6': case '7' : case '8': case '9':
        case 't': case 'f': case 'n' :
          /* And they must not be keys of the object */
          if (tokens != NULL && parser->toksuper != -1)
            {
              jsmntok_t *t = &tokens[parser->toksuper];
              if (t->type == JSMN_OBJECT || (t->type == JSMN_STRING && t->size != 0))
                return JSMN_ERROR_INVAL;
            }
#else
          /* In non-strict mode every unquoted value is a primitive */
        default:
#endif
          r = jsmn_parse_primitive(parser, js, len, tokens, num_tokens);
          if (r < 0) return r;
          count++;
          if (parser->toksuper != -1 && tokens != NULL)
            tokens[parser->toksuper].size++;
          break;

#ifdef JSMN_STRICT
          /* Unexpected char in strict mode */
        default:
          return JSMN_ERROR_INVAL;
#endif
        }
    }

  if ( parser->pos == len )
    {
          if (tokens == NULL) return JSMN_ERROR_NOMEM;
          int type = JSMN_OBJECT;
          for (i = parser->toknext - 1; i >= 0; i--)
            {
              token = &tokens[i];
              if (token->start != -1 && token->end == -1)
                {
                  if (token->type != type) return JSMN_ERROR_INVAL;
                  parser->toksuper = -1;
                  token->end = parser->pos + 1;
                  break;
                }
            }
          printf("xend %d\n", token->end);
          /* Error if unmatched closing bracket */
          if (i == -1) return JSMN_ERROR_INVAL;
          for (; i >= 0; i--)
            {
              token = &tokens[i];
              if (token->start != -1 && token->end == -1) {
                parser->toksuper = i;
                break;
              }
            }
    }

  if (tokens != NULL)
    {
      for (i = parser->toknext - 1; i >= 0; i--)
        {
          /* Unmatched opened object or array */
          if (tokens[i].start != -1 && tokens[i].end == -1) return JSMN_ERROR_PART;
        }
    }

  return count;
}

/**
 * Creates a new parser based over a given  buffer with an array of tokens
 * available.
 */
void jsmn_init(jsmn_parser *parser)
{
  parser->pos = 0;
  parser->toknext = 0;
  parser->toksuper = -1;
}

