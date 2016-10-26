#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>


typedef enum {
  NAMELIST_UNDEFINED = 0,
  NAMELIST_OBJECT    = 1,
  NAMELIST_KEY       = 2,
  NAMELIST_STRING    = 3,
  NAMELIST_WORD      = 4
} namelisttype_t;


enum namelisterr {
  NAMELIST_ERROR_INVAL = -1,   // Invalid character inside NAMELIST string/word
  NAMELIST_ERROR_PART  = -2,   // The string is not a full NAMELIST packet, more bytes expected 
  NAMELIST_ERROR_INKEY = -3,   // Invalid character inside NAMELIST key
  NAMELIST_ERROR_INTYP = -4,   // Invalid NAMELIST key type
  NAMELIST_ERROR_INOBJ = -5    // Invalid NAMELIST object
};

// NAMELIST token description.
typedef struct {
  namelisttype_t type; // type (object, key, string word)
  int start;           // start position in NAMELIST buffer
  int end;             // end position in NAMELIST buffer
} namelisttok_t;


typedef struct {
  namelisttok_t *tokens;
  unsigned int num_tokens;
  unsigned int toknext;
  unsigned int pos;
  unsigned int lineno;
} namelist_parser;

static
void namelist_init(namelist_parser *parser)
{
  parser->tokens = NULL;
  parser->num_tokens = 0;
  parser->toknext = 0;
  parser->pos = 0;
}


namelist_parser *namelist_new(void)
{
  namelist_parser *parser = (namelist_parser *) malloc(sizeof(namelist_parser));
  namelist_init(parser);

  return parser;
}


void namelist_destroy(namelist_parser *parser)
{
  if ( parser )
    {
      if ( parser->tokens ) free(parser->tokens);
      namelist_init(parser);
      free(parser);
    }
}

// Allocates a fresh unused token from the token pull.
static
namelisttok_t *namelist_alloc_token(namelist_parser *parser)
{
  const unsigned int TOK_MEM_INCR = 64;
  namelisttok_t *tok;

  if ( parser->toknext >= parser->num_tokens )
    {
      parser->num_tokens += TOK_MEM_INCR;
      parser->tokens = (namelisttok_t *) realloc(parser->tokens, sizeof(*tok) * parser->num_tokens);
      if ( parser->tokens == NULL )
        {
          fprintf(stderr, "%s: Failed to allocated more memory!", __func__);
          exit(-1);
        }
    }

  tok = &parser->tokens[parser->toknext++];
  tok->start = tok->end = -1;
  return tok;
}

// Fills token type and boundaries.
static
void namelist_fill_token(namelisttok_t *token, namelisttype_t type, int start, int end)
{
  token->type = type;
  token->start = start;
  token->end = end;
}


void namelist_new_object(namelist_parser *parser)
{
  namelisttok_t *token;
  token = namelist_alloc_token(parser);
  token->type = NAMELIST_OBJECT;
  token->start = parser->pos;
}

// Fills next available token with NAMELIST word.
static
int namelist_parse_word(namelist_parser *parser, const char *buf, size_t len)
{
  namelisttok_t *token;
  int start = parser->pos;

  for ( ; parser->pos < len && buf[parser->pos] != '\0'; parser->pos++ )
    {
      switch (buf[parser->pos])
        {
        case ':': case '=':
        case ',': case '&': case '/': 
        case '\r': case '\n': case '\t': case ' ':
          goto found;
        }

      if ( buf[parser->pos] < 32 || buf[parser->pos] >= 127 )
        {
          parser->pos = start;
          return NAMELIST_ERROR_INVAL;
        }
    }

 found:

  token = namelist_alloc_token(parser);
  namelist_fill_token(token, NAMELIST_WORD, start, parser->pos);
  parser->pos--;

  return 0;
}

// Fills next token with NAMELIST string.
static
int namelist_parse_string(namelist_parser *parser, const char *buf, size_t len)
{
  namelisttok_t *token;
  int start = parser->pos;

  parser->pos++;

  /* Skip starting quote */
  for ( ; parser->pos < len && buf[parser->pos] != '\0'; parser->pos++ )
    {
      char c = buf[parser->pos];

      /* Quote: end of string */
      if ( c == '\"' )
        {
          token = namelist_alloc_token(parser);
          namelist_fill_token(token, NAMELIST_STRING, start+1, parser->pos);
          return 0;
        }

      /* Backslash: Quoted symbol expected */
      if ( c == '\\' && parser->pos + 1 < len )
        {
          parser->pos++;
          switch (buf[parser->pos])
            {
              // Allowed escaped symbols
            case '\"': case '\\' : case 'b' :
            case 'f' : case 'r' : case 'n'  : case 't' :
              break;
              // Allows escaped symbol \uXXXX
            case 'u':
              parser->pos++;
              for ( int i = 0; i < 4 && parser->pos < len && buf[parser->pos] != '\0'; i++ )
                {
                  // If it isn't a hex character we have an error
                  if ( !((buf[parser->pos] >= 48 && buf[parser->pos] <= 57) || // 0-9
                         (buf[parser->pos] >= 65 && buf[parser->pos] <= 70) || // A-F
                         (buf[parser->pos] >= 97 && buf[parser->pos] <= 102)) ) // a-f
                    {
                      return NAMELIST_ERROR_INVAL;
                    }
                  parser->pos++;
                }
              parser->pos--;
              break;
              // Unexpected symbol
            default:
              return NAMELIST_ERROR_INVAL;
            }
        }
    }
  
  parser->pos = start;
  return NAMELIST_ERROR_PART;

  return 0;
}

static
int namelist_check_keyname(const char *buf, namelisttok_t *t)
{
  switch ((int)(t->type))
    {
    case (int)NAMELIST_STRING:
      for ( int i = t->start; i < t->end; ++i )
        {
          switch(buf[i])
            {
            case '\t': case ' ': return NAMELIST_ERROR_INKEY;
            }
        }
    case (int)NAMELIST_WORD:
      t->type = NAMELIST_KEY;
      break;
    default:
      return NAMELIST_ERROR_INTYP;
      break;
    }
  
  return 0;
}


int namelist_parse(namelist_parser *parser, const char *buf, size_t len)
{
  int status = 0;
  namelisttok_t *token;

  parser->lineno = 1;

  for ( ; parser->pos < len && buf[parser->pos] != '\0'; parser->pos++ )
    {
      char c = buf[parser->pos];
      switch (c)
        {
        case '&':
          namelist_new_object(parser);
          break;
        case '/':
          for ( int i = parser->toknext - 1; i >= 0; i-- )
            {
              token = &parser->tokens[i];
              if ( token->start != -1 && token->end == -1 )
                {
                  if ( token->type != NAMELIST_OBJECT ) return NAMELIST_ERROR_INOBJ;
                  token->end = parser->pos + 1;
                  break;
                }
            }
          break;
        case '\t': case ' ':
          break;
        case '\r':
          parser->lineno++;
          if ( parser->pos+1 < len && buf[parser->pos+1] == '\n' ) parser->pos++;
          break;
        case '\n':
          parser->lineno++;
          break;
        case ',':
          break;
        case '#': case '!': // Skip to end of line
          for ( ; parser->pos < len && buf[parser->pos] != '\0'; parser->pos++ )
            if ( buf[parser->pos] == '\r' || buf[parser->pos] == '\n' )
              {
                parser->pos--;
                break;
              }
          break;
        case ':': case '=':
          status = namelist_check_keyname(buf, &parser->tokens[parser->toknext-1]);
          break;
        case '\"':
          status = namelist_parse_string(parser, buf, len);
          break;
        default:
          status = namelist_parse_word(parser, buf, len);
          break;
        }

      if ( status ) return status;
    }

  return status;
}


void namelist_dump(namelist_parser *parser, const char *buf)
{
  unsigned int ntok = parser->toknext;
  printf("Number of tokens %d\n", ntok);

  for ( unsigned int it = 0; it < ntok; ++it )
    {
      namelisttok_t *t = &parser->tokens[it];
      int length = t->end - t->start;
      const char *start = buf+t->start;
      printf("Token %u", it+1);
      if ( t->type == NAMELIST_OBJECT )
        {
          printf(" NAMELIST=");
          if ( length > 80 ) length = 80;
          printf("'%.*s'", length, start);
        }
      else if ( t->type == NAMELIST_KEY )
        {
          printf(" KEY=");
          printf("'%.*s'", length, start);
        }
      else if ( t->type == NAMELIST_WORD )
        {
          printf(" WORD=");
          printf("'%.*s'", length, start);
        }
      else if ( t->type == NAMELIST_STRING )
        {
          printf(" STRING=");
          printf("'%.*s'", length, start);
        }
      printf("\n");
    }
}

int main(int argc, char *argv[])
{
  if ( argc != 2 )
    {
      fprintf(stderr, "Usage: %s namelist\n", argv[0]);
      return -1;
    }

  const char *filename = argv[1];
  printf("Parse namelist %s:\n", filename);

  struct stat sbuf;
  size_t filesize = (stat(filename, &sbuf) == 0) ? sbuf.st_size : 0;

  if ( filesize == 0 )
    {
      fprintf(stderr, "Empty table file: %s\n", filename);
      return -1;
    }

  FILE *fp = fopen(filename, "r");
  if ( fp == NULL )
    {
      fprintf(stderr, "Open failed on %s: %s\n", filename, strerror(errno));
      return -1;
    }

  char *buffer = (char*) malloc(filesize);
  size_t nitems = fread(buffer, 1, filesize, fp);

  fclose(fp);

  if ( nitems != filesize )
    {
      fprintf(stderr, "Read failed on %s!\n", filename);
      return -1;
    }

  namelist_parser *p = namelist_new();

  int status = namelist_parse(p, buffer, filesize);
  printf("Processed number of lines: %d\n", p->lineno-1);
  if ( status == NAMELIST_ERROR_INVAL ) fprintf(stderr, "Namelist parser error: invalid character in %s (line=%d character='%c')!\n", filename, p->lineno, buffer[p->pos]);
  if ( status == NAMELIST_ERROR_PART  ) fprintf(stderr, "Namelist parser error: end of string not found in %s (line=%d)!\n", filename, p->lineno);
  if ( status == NAMELIST_ERROR_INKEY ) fprintf(stderr, "Namelist parser error: invalid key word in %s (line=%d)!\n", filename, p->lineno);
  if ( status == NAMELIST_ERROR_INTYP ) fprintf(stderr, "Namelist parser error: invalid key word type in %s (line=%d)!\n", filename, p->lineno);
  if ( status == NAMELIST_ERROR_INOBJ ) fprintf(stderr, "Namelist parser error: invalid object in %s (line=%d)!\n", filename, p->lineno);
  namelist_dump(p, buffer);

  namelist_destroy(p);

  free(buffer);

  return 0;
}
