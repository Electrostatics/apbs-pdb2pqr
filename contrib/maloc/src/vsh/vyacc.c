/* A Bison parser, made from ../../../src/vsh/vyacc.y, by GNU bison 1.75.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Written by Richard Stallman by simplifying the original so called
   ``semantic'' parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON	1

/* Pure parsers.  */
#define YYPURE	0

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     IF = 258,
     THEN = 259,
     ELSE = 260,
     ELIF = 261,
     FI = 262,
     CASE = 263,
     IN = 264,
     ESAC = 265,
     FOR = 266,
     WHILE = 267,
     UNTIL = 268,
     DO = 269,
     DONE = 270,
     SEMI_SEMI = 271,
     AND_AND = 272,
     OR_OR = 273,
     LESS_LESS = 274,
     GREATER_GREATER = 275,
     LESS_AND = 276,
     GREATER_AND = 277,
     AND_GREATER = 278,
     LESS_GREATER = 279,
     GREATER_BAR = 280,
     LESS_LESS_MINUS = 281,
     WORD = 282,
     ASSIGNMENT_WORD = 283,
     SELECT = 284,
     FUNCTION = 285,
     BANG = 286,
     vsh_EOF = 287
   };
#endif
#define IF 258
#define THEN 259
#define ELSE 260
#define ELIF 261
#define FI 262
#define CASE 263
#define IN 264
#define ESAC 265
#define FOR 266
#define WHILE 267
#define UNTIL 268
#define DO 269
#define DONE 270
#define SEMI_SEMI 271
#define AND_AND 272
#define OR_OR 273
#define LESS_LESS 274
#define GREATER_GREATER 275
#define LESS_AND 276
#define GREATER_AND 277
#define AND_GREATER 278
#define LESS_GREATER 279
#define GREATER_BAR 280
#define LESS_LESS_MINUS 281
#define WORD 282
#define ASSIGNMENT_WORD 283
#define SELECT 284
#define FUNCTION 285
#define BANG 286
#define vsh_EOF 287




/* Copy the first part of user declarations.  */


#include "vsh_p.h"

VEMBED(rcsid="$Id: vyacc.c,v 1.10 2003/07/05 20:44:12 mholst Exp $")

VPUBLIC void yyerror(const char *errmsg)
{
    fprintf(stderr, "Vsh: %s at '%s'\n", errmsg, yytext);
    if (!Vsh_isInteractive(Vsh_thee)) {
        exit(1);
    }
}

VPUBLIC void report_error(const char *format, ...)
{
    fprintf(stderr,"report_error: ERROR:\n");
}

VPUBLIC void programming_error(const char *format, ...)
{
    fprintf(stderr,"programming_error: ERROR:\n");
}

VPUBLIC int yywrap(void)
{
    if (Vsh_thee->cinUnit == stdin) {
        Vnm_print(1,"%s",Vsh_thee->PR_EXIT);
        Vnm_print(1,"%s",VNEWLINE_STRING);
    }
    cmdKey = 2;
    return 1;
}

/*
 * Reverse the chain of structures in LIST.  Output the new head
 * of the chain.  You should always assign the output value of this
 * function to something, or you will lose the chain.
 */
GENERIC_LIST *reverse_list (GENERIC_LIST *list)
{
    register GENERIC_LIST *next, *prev = (GENERIC_LIST *)NULL;

    while (list) {
        next = list->next;
        list->next = prev;
        prev = list;
        list = next;
    }
    return (prev);
}

WORD_LIST* make_word_list(WORD_DESC *word, WORD_LIST *link)
{
    WORD_LIST *temp;

    temp = (WORD_LIST *)malloc (sizeof (WORD_LIST));
    temp->word = word;
    temp->next = link;
    return temp;
}

COMMAND *make_bare_simple_command(void)
{
    COMMAND *command;
    SIMPLE_COM *temp = (SIMPLE_COM *)malloc (sizeof (SIMPLE_COM));

    temp->flags = 0;
    /* temp->line = line_number; */
    temp->line = 0;
    temp->words = (WORD_LIST *)NULL;
    temp->redirects = (REDIRECT *)NULL;
    command = (COMMAND *)malloc (sizeof (COMMAND));
    command->type = cm_simple;
    command->redirects = (REDIRECT *)NULL;
    command->flags = 0;
    command->value.Simple = temp;
    return command;
}

/*
 * Generate a REDIRECT from SOURCE, DEST, and INSTRUCTION. 
 * INSTRUCTION is the instruction type, SOURCE is a file descriptor,
 * and DEST is a file descriptor or a WORD_DESC*.
 */
REDIRECT *make_redirection (int source,
    enum r_instruction instruction, REDIRECTEE dest_and_filename)
{
    REDIRECT *temp = (REDIRECT *)malloc (sizeof (REDIRECT));

    /* First do the common cases. */
    temp->redirector = source;
    temp->redirectee = dest_and_filename;
    temp->instruction = instruction;
    temp->flags = 0;
    temp->next = (REDIRECT *)NULL;

    switch (instruction) {

      case r_output_direction:	/* >foo */
      case r_output_force:	/* >| foo */
        temp->flags = O_TRUNC | O_WRONLY | O_CREAT;
        break;

      case r_input_direction:	/* <foo */
      case r_inputa_direction:	/* foo & makes this. */
        temp->flags = O_RDONLY;
        break;

      case r_appending_to:	/* >>foo */
        temp->flags = O_APPEND | O_WRONLY | O_CREAT;
        break;

      case r_deblank_reading_until: /* <<-foo */
      case r_reading_until:	/* << foo */
        break;

      case r_duplicating_input:		/* 1<&2 */
      case r_duplicating_output:		/* 1>&2 */
      case r_close_this:			/* <&- */
      case r_duplicating_input_word:	/* 1<&$foo */
      case r_duplicating_output_word:	/* 1>&$foo */
        break;
    
      case r_err_and_out:		/* command &>filename */
        temp->flags = O_TRUNC | O_WRONLY | O_CREAT;
        break;

      case r_input_output:
        temp->flags = O_RDWR | O_CREAT;
        break;

      default:
        fprintf(stderr,"We seem to have a problem...\n");
        abort ();
        break;
      }
    return (temp);
}

/*
 * Return a command which is the connection of the word or redirection
 * in ELEMENT, and the command * or NULL in COMMAND.
 */
COMMAND *make_simple_command (ELEMENT element, COMMAND *command) {
    /* If we are starting from scratch, then make the initial command
       structure.  Also note that we have to fill in all the slots, since
       malloc doesn't return zeroed space. */
    if (!command) command = make_bare_simple_command ();

    if (element.word) {
        WORD_LIST *tw = (WORD_LIST *)malloc (sizeof (WORD_LIST));
        tw->word = element.word;
        tw->next = command->value.Simple->words;
        command->value.Simple->words = tw;
    } else {
        REDIRECT *r = element.redirect;
        /* Due to the way <> is implemented, there may be more than a single
           redirection in element.redirect.  We just follow the chain as far
           as it goes, and hook onto the end. */
        while (r->next) r = r->next;
        r->next = command->value.Simple->redirects;
        command->value.Simple->redirects = element.redirect;
    }
    return (command);
}

COMMAND* make_command (enum command_type type, SIMPLE_COM *pointer)
{
    COMMAND *temp;

    temp = (COMMAND *)malloc (sizeof (COMMAND));
    temp->type = type;
    temp->value.Simple = pointer;
    temp->value.Simple->flags = 0;
    temp->flags = 0;
    temp->redirects = (REDIRECT *)NULL;
    return temp;
}

void dispose_command (COMMAND *command);
void dispose_word (WORD_DESC *word);
void dispose_words (WORD_LIST *list);
void dispose_word_array (char **array);
void dispose_redirects (REDIRECT *list);

/* Dispose of the command structure passed. */
void dispose_command (COMMAND *command)
{
  if (!command) return;

  if (command->redirects)
    dispose_redirects (command->redirects);

  switch (command->type)
    {
    case cm_for:
      {
	register FOR_COM *c = command->value.For;
	dispose_word (c->name);
	dispose_words (c->map_list);
	dispose_command (c->action);
	free (c);
	break;
      }

    case cm_group:
      {
	dispose_command (command->value.Group->command);
	free (command->value.Group);
	break;
      }

    case cm_case:
      {
	register CASE_COM *c = command->value.Case;
	PATTERN_LIST *t, *p = c->clauses;

	dispose_word (c->word);

	while (p)
	  {
	    dispose_words (p->patterns);
	    dispose_command (p->action);
	    t = p;
	    p = p->next;
	    free (t);
	  }
	free (c);
	break;
      }

    case cm_until:
    case cm_while:
      {
	register WHILE_COM *c = command->value.While;

	dispose_command (c->test);
	dispose_command (c->action);
	free (c);
	break;
      }

    case cm_if:
      {
	register IF_COM *c = command->value.If;
	dispose_command (c->test);
	dispose_command (c->true_case);
	dispose_command (c->false_case);
	free (c);
	break;
      }

    case cm_simple:
      {
	register SIMPLE_COM *c = command->value.Simple;
	dispose_words (c->words);
	dispose_redirects (c->redirects);
	free (c);
	break;
      }

    case cm_connection:
      {
	register CONNECTION *c = command->value.Connection;
	dispose_command (c->first);
	dispose_command (c->second);
	free (c);
	break;
      }

    case cm_function_def:
      {
	register FUNCTION_DEF *c = command->value.Function_def;
	dispose_word (c->name);
	dispose_command (c->command);
	free (c);
	break;
      }

    default:
      report_error ("Attempt free unknown command type `%d'.\n", command->type);
      break;
    }
  free (command);
}

/* How to free a WORD_DESC. */
void dispose_word (WORD_DESC *word)
{
  if (word->word)
    free (word->word);
  free (word);
}

/* How to get rid of a linked list of words.  A WORD_LIST. */
void dispose_words (WORD_LIST *list)
{
  WORD_LIST *t;
  while (list)
    {
      t = list;
      list = list->next;
      dispose_word (t->word);
      free (t);
    }
}

/* How to dispose of an array of pointers to char. */
void dispose_word_array (char **array)
{
  register int count;

  for (count = 0; array[count]; count++)
    free (array[count]);

  free (array);
}

/* How to dispose of an list of redirections.  A REDIRECT. */
void dispose_redirects (REDIRECT *list)
{
  register REDIRECT *t;

  while (list)
    {
      t = list;
      list = list->next;
      switch (t->instruction)
	{
	case r_reading_until:
	case r_deblank_reading_until:
	  free (t->here_doc_eof);
	  /* ... */
	case r_output_direction:
	case r_input_direction:
	case r_inputa_direction:
	case r_appending_to:
	case r_err_and_out:
	case r_input_output:
	case r_output_force:
	case r_duplicating_input_word:
	case r_duplicating_output_word:
	  dispose_word (t->redirectee.filename);
	  break;
	}
      free (t);
    }
}

/* Reverse the word list and redirection list in the simple command
   has just been parsed.  It seems simpler to do this here the one
   time then by any other method that I can think of. */
COMMAND *clean_simple_command (COMMAND *command)
{
  if (command->type != cm_simple)
    {
      programming_error
    ("clean_simple_command () got a command with type %d.", command->type);
    }
  else
    {
      command->value.Simple->words =
    REVERSE_LIST (command->value.Simple->words, WORD_LIST *);
      command->value.Simple->redirects =
    REVERSE_LIST (command->value.Simple->redirects, REDIRECT *);
    }

  return (command);
}

static REDIRECTEE redir;

static WORD_LIST *wordTmp;
static WORD_DESC *wTmp;
static char buf[VMAX_BUFSIZE];



/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

#ifndef YYSTYPE
#line 413 "../../../src/vsh/vyacc.y"
typedef union {
    WORD_DESC *word;         /* the word that we read. */
    /* int number; */        /* number we saw */
    WORD_LIST *word_list;    /* a sequence of white-space separated words */
    COMMAND *command;        /* a complete command */
    REDIRECT *redirect;      /* redirect i/o info */
    ELEMENT element;         /* base element */
    PATTERN_LIST *pattern;   /* a case pattern */
} yystype;
/* Line 193 of /usr/share/bison/yacc.c.  */
#line 525 "y.tab.c"
# define YYSTYPE yystype
# define YYSTYPE_IS_TRIVIAL 1
#endif

#ifndef YYLTYPE
typedef struct yyltype
{
  int first_line;
  int first_column;
  int last_line;
  int last_column;
} yyltype;
# define YYLTYPE yyltype
# define YYLTYPE_IS_TRIVIAL 1
#endif

/* Copy the second part of user declarations.  */


/* Line 213 of /usr/share/bison/yacc.c.  */
#line 546 "y.tab.c"

#if ! defined (yyoverflow) || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# if YYSTACK_USE_ALLOCA
#  define YYSTACK_ALLOC alloca
# else
#  ifndef YYSTACK_USE_ALLOCA
#   if defined (alloca) || defined (_ALLOCA_H)
#    define YYSTACK_ALLOC alloca
#   else
#    ifdef __GNUC__
#     define YYSTACK_ALLOC __builtin_alloca
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
# else
#  if defined (__STDC__) || defined (__cplusplus)
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   define YYSIZE_T size_t
#  endif
#  define YYSTACK_ALLOC malloc
#  define YYSTACK_FREE free
# endif
#endif /* ! defined (yyoverflow) || YYERROR_VERBOSE */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
	 || (YYLTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAX (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short) + sizeof (YYSTYPE))				\
      + YYSTACK_GAP_MAX)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  register YYSIZE_T yyi;		\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];	\
	}					\
      while (0)
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAX;	\
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (0)

#endif

#if defined (__STDC__) || defined (__cplusplus)
   typedef signed char yysigned_char;
#else
   typedef short yysigned_char;
#endif

/* YYFINAL -- State number of the termination state. */
#define YYFINAL  47
#define YYLAST   423

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  43
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  32
/* YYNRULES -- Number of rules. */
#define YYNRULES  91
/* YYNRULES -- Number of states. */
#define YYNSTATES  181

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   287

#define YYTRANSLATE(X) \
  ((unsigned)(X) <= YYMAXUTOK ? yytranslate[X] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const unsigned char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      34,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    37,     2,     2,     2,     2,    32,     2,
      39,    40,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,    33,
       2,     2,    38,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    41,    36,    42,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    35
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const unsigned short yyprhs[] =
{
       0,     0,     3,     6,     8,    10,    13,    15,    18,    21,
      26,    31,    35,    39,    41,    44,    49,    51,    54,    56,
      60,    64,    68,    73,    78,    83,    88,    93,    95,    98,
     101,   104,   106,   109,   111,   113,   115,   118,   120,   123,
     125,   127,   129,   131,   133,   135,   137,   139,   141,   143,
     145,   149,   153,   159,   166,   172,   180,   187,   192,   199,
     205,   212,   220,   227,   229,   232,   237,   242,   248,   254,
     256,   259,   265,   271,   278,   285,   287,   291,   298,   305,
     313,   321,   332,   343,   349,   355,   356,   359,   361,   363,
     365,   366
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const yysigned_char yyrhs[] =
{
      44,     0,    -1,    45,    34,    -1,    34,    -1,    35,    -1,
       1,    74,    -1,    46,    -1,    46,    32,    -1,    46,    33,
      -1,    46,    17,    74,    46,    -1,    46,    18,    74,    46,
      -1,    46,    32,    46,    -1,    46,    33,    46,    -1,    47,
      -1,    37,    47,    -1,    47,    36,    74,    47,    -1,    53,
      -1,    74,    49,    -1,    50,    -1,    50,    34,    74,    -1,
      50,    32,    74,    -1,    50,    33,    74,    -1,    50,    17,
      74,    50,    -1,    50,    18,    74,    50,    -1,    50,    32,
      74,    50,    -1,    50,    33,    74,    50,    -1,    50,    34,
      74,    50,    -1,    47,    -1,    37,    47,    -1,    38,    27,
      -1,    20,    27,    -1,    51,    -1,    52,    51,    -1,    54,
      -1,    55,    -1,    57,    -1,    54,    57,    -1,    56,    -1,
      56,    52,    -1,    69,    -1,    63,    -1,    70,    -1,    71,
      -1,    61,    -1,    58,    -1,    59,    -1,    60,    -1,    27,
      -1,    28,    -1,    51,    -1,    39,    48,    40,    -1,    41,
      48,    42,    -1,    27,    39,    40,    74,    59,    -1,    27,
      39,    40,    74,    59,    52,    -1,     3,    48,     4,    48,
       7,    -1,     3,    48,     4,    48,     5,    48,     7,    -1,
       3,    48,     4,    48,    62,     7,    -1,     6,    48,     4,
      48,    -1,     6,    48,     4,    48,     5,    48,    -1,     6,
      48,     4,    48,    62,    -1,     8,    27,    74,     9,    74,
      10,    -1,     8,    27,    74,     9,    66,    74,    10,    -1,
       8,    27,    74,     9,    64,    10,    -1,    65,    -1,    66,
      65,    -1,    74,    68,    40,    48,    -1,    74,    68,    40,
      74,    -1,    74,    39,    68,    40,    48,    -1,    74,    39,
      68,    40,    74,    -1,    67,    -1,    66,    67,    -1,    74,
      68,    40,    48,    16,    -1,    74,    68,    40,    74,    16,
      -1,    74,    39,    68,    40,    48,    16,    -1,    74,    39,
      68,    40,    74,    16,    -1,    27,    -1,    68,    36,    27,
      -1,    11,    27,    74,    14,    48,    15,    -1,    11,    27,
      74,    41,    48,    42,    -1,    11,    27,    33,    74,    14,
      48,    15,    -1,    11,    27,    33,    74,    41,    48,    42,
      -1,    11,    27,    74,     9,    72,    73,    74,    14,    48,
      15,    -1,    11,    27,    74,     9,    72,    73,    74,    41,
      48,    42,    -1,    12,    48,    14,    48,    15,    -1,    13,
      48,    14,    48,    15,    -1,    -1,    72,    27,    -1,    34,
      -1,    33,    -1,    35,    -1,    -1,    74,    34,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned short yyrline[] =
{
       0,   471,   471,   496,   501,   507,   526,   530,   534,   540,
     544,   548,   552,   556,   560,   566,   570,   582,   588,   592,
     596,   600,   606,   610,   614,   618,   622,   626,   630,   642,
     647,   687,   691,   703,   707,   713,   717,   723,   727,   733,
     737,   741,   745,   749,   753,   757,   761,   767,   772,   777,
     790,   796,   802,   806,   818,   822,   826,   832,   836,   840,
     852,   856,   860,   866,   870,   876,   880,   884,   888,   894,
     898,   904,   908,   912,   916,   922,   926,   938,   942,   946,
     950,   954,   958,   970,   982,   994,   998,  1004,  1007,  1010,
    1015,  1018
};
#endif

#if YYDEBUG || YYERROR_VERBOSE
/* YYTNME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals. */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "IF", "THEN", "ELSE", "ELIF", "FI", "CASE", 
  "IN", "ESAC", "FOR", "WHILE", "UNTIL", "DO", "DONE", "SEMI_SEMI", 
  "AND_AND", "OR_OR", "LESS_LESS", "GREATER_GREATER", "LESS_AND", 
  "GREATER_AND", "AND_GREATER", "LESS_GREATER", "GREATER_BAR", 
  "LESS_LESS_MINUS", "WORD", "ASSIGNMENT_WORD", "SELECT", "FUNCTION", 
  "BANG", "'&'", "';'", "'\\n'", "vsh_EOF", "'|'", "'!'", "'>'", "'('", 
  "')'", "'{'", "'}'", "$accept", "inputunit", "simple_list", 
  "simple_list1", "pipeline", "list", "list0", "list1", "redirection", 
  "redirections", "command", "simple_command", "shell_command", 
  "shell_command_1", "simple_command_element", "subshell", 
  "group_command", "function_def", "if_command", "elif_clause", 
  "case_command", "case_clause_1", "pattern_list_1", 
  "case_clause_sequence", "pattern_list", "pattern", "for_command", 
  "while_command", "until_command", "words", "term", "newlines", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const unsigned short yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,    38,    59,    10,   287,   124,    33,    62,    40,
      41,   123,   125
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
       0,    43,    44,    44,    44,    44,    45,    45,    45,    46,
      46,    46,    46,    46,    46,    47,    47,    48,    49,    49,
      49,    49,    50,    50,    50,    50,    50,    50,    50,    51,
      51,    52,    52,    53,    53,    54,    54,    55,    55,    56,
      56,    56,    56,    56,    56,    56,    56,    57,    57,    57,
      58,    59,    60,    60,    61,    61,    61,    62,    62,    62,
      63,    63,    63,    64,    64,    65,    65,    65,    65,    66,
      66,    67,    67,    67,    67,    68,    68,    69,    69,    69,
      69,    69,    69,    70,    71,    72,    72,    73,    73,    73,
      74,    74
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const unsigned char yyr2[] =
{
       0,     2,     2,     1,     1,     2,     1,     2,     2,     4,
       4,     3,     3,     1,     2,     4,     1,     2,     1,     3,
       3,     3,     4,     4,     4,     4,     4,     1,     2,     2,
       2,     1,     2,     1,     1,     1,     2,     1,     2,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       3,     3,     5,     6,     5,     7,     6,     4,     6,     5,
       6,     7,     6,     1,     2,     4,     4,     5,     5,     1,
       2,     5,     5,     6,     6,     1,     3,     6,     6,     7,
       7,    10,    10,     5,     5,     0,     2,     1,     1,     1,
       0,     2
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned char yydefact[] =
{
       0,    90,    90,     0,     0,    90,    90,     0,    47,    48,
       3,     4,     0,     0,    90,    90,     0,     0,     6,    13,
      49,    16,    33,    34,    37,    35,    44,    45,    46,    43,
      40,    39,    41,    42,     5,     0,     0,    90,    90,     0,
       0,    30,     0,    14,    29,     0,     0,     1,     2,    90,
      90,     7,     8,    90,    47,    36,    31,    38,    91,    90,
       0,    27,    17,    18,     0,    90,     0,    90,    90,    90,
      50,    51,     0,     0,    11,    12,     0,    32,     0,    28,
      90,    90,    90,    90,    90,    90,     0,    85,    90,    90,
       0,     0,     0,     9,    10,     0,     0,    15,    90,    90,
      54,     0,     0,     0,    20,    21,    19,     0,    63,    90,
      69,     0,    90,    90,     0,     0,     0,    83,    84,    52,
       0,     0,    56,    22,    23,    24,    25,    26,    62,    64,
      70,     0,    60,    75,     0,     0,     0,     0,    86,    88,
      87,    89,    90,    77,    78,    53,    55,    90,    90,    90,
      90,    61,     0,     0,    90,    79,    80,     0,    57,     0,
       0,     0,    90,    76,    65,    66,    90,    90,    90,    59,
      67,    68,    71,    72,     0,     0,    58,    73,    74,    81,
      82
};

/* YYDEFGOTO[NTERM-NUM]. */
static const short yydefgoto[] =
{
      -1,    16,    17,    74,    61,    35,    62,    63,    20,    57,
      21,    22,    23,    24,    25,    26,    27,    28,    29,   101,
      30,   107,   108,   109,   110,   135,    31,    32,    33,   114,
     142,    36
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -92
static const short yypact[] =
{
     163,   -92,   -92,     1,    13,   -92,   -92,    18,    11,   -92,
     -92,   -92,   382,    26,   -92,   -92,    17,    33,    78,    43,
     -92,   -92,     9,   -92,    14,   -92,   -92,   -92,   -92,   -92,
     -92,   -92,   -92,   -92,    42,    88,   264,   -92,    52,    68,
      91,   -92,    76,    43,   -92,    77,    82,   -92,   -92,   -92,
     -92,   328,   328,   -92,   -92,   -92,   -92,    14,   -92,   -92,
     382,    43,   -92,    40,    -4,   -92,     7,   -92,   -92,   -92,
     -92,   -92,   296,   296,    60,    60,   350,   -92,    59,    43,
     -92,   -92,   -92,   -92,   -92,   -92,   -10,   -92,   -92,   -92,
     103,   104,   -16,   -92,   -92,   328,   328,    43,   -92,   -92,
     -92,   113,   264,   264,   264,   264,   264,   116,   -92,   -92,
     -92,    -7,   -92,   -92,    79,   112,    86,   -92,   -92,    14,
     122,   126,   -92,   -92,   -92,    81,    81,    81,   -92,   -92,
     -92,    12,   -92,   -92,   105,   -17,   118,    89,   -92,   -92,
     -92,   -92,   -92,   -92,   -92,    14,   -92,   -92,   -92,   -92,
     -92,   -92,     2,   107,   -92,   -92,   -92,    -8,    97,   264,
     264,   264,   -92,   -92,   119,   200,   -92,   -92,   -92,   -92,
     120,   232,   -92,   -92,   123,    95,   -92,   -92,   -92,   -92,
     -92
};

/* YYPGOTO[NTERM-NUM].  */
static const yysigned_char yypgoto[] =
{
     -92,   -92,   -92,     8,    49,    -5,   -92,   -91,   -22,    20,
     -92,   -92,   -92,   -92,   121,   -92,    48,   -92,   -92,   -12,
     -92,   -92,    32,   -92,    38,    16,   -92,   -92,   -92,   -92,
     -92,     6
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, parse error.  */
#define YYTABLE_NINF -1
static const unsigned char yytable[] =
{
      39,    40,    56,   132,   112,    85,   166,    34,    18,    45,
      46,   123,   124,   125,   126,   127,    87,    47,    58,   153,
     133,    88,   151,   154,    58,    15,    58,    58,    37,     7,
      58,   113,   134,   167,     7,    77,    54,     9,   153,   133,
      38,    58,   162,    64,    66,    41,    58,    13,    89,    19,
      42,   134,    13,    44,    78,    72,    73,    80,    81,    76,
      75,    43,    90,    91,    98,    99,   100,    48,   125,   126,
     127,    86,    82,    83,    84,    92,    58,    49,    50,    53,
      93,    94,    67,   115,   116,    65,   102,   103,   104,   105,
     106,   111,    59,   120,   121,    49,    50,    56,    80,    81,
      19,    19,   168,    99,    75,    68,   138,   136,   137,    79,
      51,    52,   139,   140,   141,   131,    69,    70,   117,   118,
     122,    19,    19,    77,    71,    97,   128,   143,   144,   146,
     147,   156,   133,   155,   163,   172,   177,   180,   179,   145,
     119,   129,   158,    55,    19,    19,   169,   130,   157,   164,
     152,     0,     0,     0,   159,   160,   161,   170,     0,     0,
     165,   174,   175,   176,     1,     0,     2,     0,   171,     0,
       0,     3,     0,     0,     4,     5,     6,     0,     0,     0,
       0,     0,     0,     7,     0,     0,     0,     0,     0,     0,
       8,     9,     0,     0,     0,     0,     0,    10,    11,     0,
      12,    13,    14,     2,    15,     0,     0,     0,     3,     0,
       0,     4,     5,     6,     0,     0,   173,     0,     0,     0,
       7,     0,     0,     0,     0,     0,     0,     8,     9,     0,
       0,     0,     0,     0,    58,     2,     0,    60,    13,    14,
       3,    15,     0,     4,     5,     6,     0,     0,   178,     0,
       0,     0,     7,     0,     0,     0,     0,     0,     0,     8,
       9,     0,     0,     0,     0,     0,    58,     2,     0,    60,
      13,    14,     3,    15,     0,     4,     5,     6,     0,     0,
       0,     0,     0,     0,     7,     0,     0,     0,     0,     0,
       0,     8,     9,     0,     0,     0,     0,     0,    58,     2,
       0,    60,    13,    14,     3,    15,     0,     4,     5,     6,
       0,     0,     0,     0,     0,     0,     7,     0,     0,     0,
       0,     0,     0,     8,     9,     0,     0,     0,     0,     0,
      58,     2,     0,    12,    13,    14,     3,    15,     0,     4,
       5,     6,     0,     0,     0,     0,     0,     0,     7,     0,
       0,     0,     0,     2,     0,     8,     9,     0,     3,     0,
       0,     4,     5,     6,     0,    12,    13,    14,     0,    15,
       7,     0,     0,     0,     0,     0,     0,     8,     9,     0,
       0,     0,     0,     0,    58,     2,     0,     0,    13,    14,
       3,    15,     0,     4,     5,     6,     0,     0,     0,     0,
       0,     0,     7,     0,     0,     0,     0,     0,     0,     8,
       9,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      13,    14,     0,    15
};

static const short yycheck[] =
{
       5,     6,    24,    10,    14,     9,    14,     1,     0,    14,
      15,   102,   103,   104,   105,   106,     9,     0,    34,    36,
      27,    14,    10,    40,    34,    41,    34,    34,    27,    20,
      34,    41,    39,    41,    20,    57,    27,    28,    36,    27,
      27,    34,    40,    37,    38,    27,    34,    38,    41,     0,
      39,    39,    38,    27,    59,    49,    50,    17,    18,    53,
      52,    12,    67,    68,     5,     6,     7,    34,   159,   160,
     161,    65,    32,    33,    34,    69,    34,    17,    18,    36,
      72,    73,    14,    88,    89,    33,    80,    81,    82,    83,
      84,    85,     4,    98,    99,    17,    18,   119,    17,    18,
      51,    52,     5,     6,    96,    14,    27,   112,   113,    60,
      32,    33,    33,    34,    35,   109,    40,    40,    15,    15,
       7,    72,    73,   145,    42,    76,    10,    15,    42,     7,
       4,    42,    27,    15,    27,    16,    16,    42,    15,   119,
      92,   109,   147,    22,    95,    96,   158,   109,   142,   154,
     134,    -1,    -1,    -1,   148,   149,   150,   162,    -1,    -1,
     154,   166,   167,   168,     1,    -1,     3,    -1,   162,    -1,
      -1,     8,    -1,    -1,    11,    12,    13,    -1,    -1,    -1,
      -1,    -1,    -1,    20,    -1,    -1,    -1,    -1,    -1,    -1,
      27,    28,    -1,    -1,    -1,    -1,    -1,    34,    35,    -1,
      37,    38,    39,     3,    41,    -1,    -1,    -1,     8,    -1,
      -1,    11,    12,    13,    -1,    -1,    16,    -1,    -1,    -1,
      20,    -1,    -1,    -1,    -1,    -1,    -1,    27,    28,    -1,
      -1,    -1,    -1,    -1,    34,     3,    -1,    37,    38,    39,
       8,    41,    -1,    11,    12,    13,    -1,    -1,    16,    -1,
      -1,    -1,    20,    -1,    -1,    -1,    -1,    -1,    -1,    27,
      28,    -1,    -1,    -1,    -1,    -1,    34,     3,    -1,    37,
      38,    39,     8,    41,    -1,    11,    12,    13,    -1,    -1,
      -1,    -1,    -1,    -1,    20,    -1,    -1,    -1,    -1,    -1,
      -1,    27,    28,    -1,    -1,    -1,    -1,    -1,    34,     3,
      -1,    37,    38,    39,     8,    41,    -1,    11,    12,    13,
      -1,    -1,    -1,    -1,    -1,    -1,    20,    -1,    -1,    -1,
      -1,    -1,    -1,    27,    28,    -1,    -1,    -1,    -1,    -1,
      34,     3,    -1,    37,    38,    39,     8,    41,    -1,    11,
      12,    13,    -1,    -1,    -1,    -1,    -1,    -1,    20,    -1,
      -1,    -1,    -1,     3,    -1,    27,    28,    -1,     8,    -1,
      -1,    11,    12,    13,    -1,    37,    38,    39,    -1,    41,
      20,    -1,    -1,    -1,    -1,    -1,    -1,    27,    28,    -1,
      -1,    -1,    -1,    -1,    34,     3,    -1,    -1,    38,    39,
       8,    41,    -1,    11,    12,    13,    -1,    -1,    -1,    -1,
      -1,    -1,    20,    -1,    -1,    -1,    -1,    -1,    -1,    27,
      28,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      38,    39,    -1,    41
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
       0,     1,     3,     8,    11,    12,    13,    20,    27,    28,
      34,    35,    37,    38,    39,    41,    44,    45,    46,    47,
      51,    53,    54,    55,    56,    57,    58,    59,    60,    61,
      63,    69,    70,    71,    74,    48,    74,    27,    27,    48,
      48,    27,    39,    47,    27,    48,    48,     0,    34,    17,
      18,    32,    33,    36,    27,    57,    51,    52,    34,     4,
      37,    47,    49,    50,    74,    33,    74,    14,    14,    40,
      40,    42,    74,    74,    46,    46,    74,    51,    48,    47,
      17,    18,    32,    33,    34,     9,    74,     9,    14,    41,
      48,    48,    74,    46,    46,    32,    33,    47,     5,     6,
       7,    62,    74,    74,    74,    74,    74,    64,    65,    66,
      67,    74,    14,    41,    72,    48,    48,    15,    15,    59,
      48,    48,     7,    50,    50,    50,    50,    50,    10,    65,
      67,    74,    10,    27,    39,    68,    48,    48,    27,    33,
      34,    35,    73,    15,    42,    52,     7,     4,    32,    33,
      34,    10,    68,    36,    40,    15,    42,    74,    48,    74,
      74,    74,    40,    27,    48,    74,    14,    41,     5,    62,
      48,    74,    16,    16,    48,    48,    48,    16,    16,    15,
      42
};

#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# if defined (__STDC__) || defined (__cplusplus)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# endif
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		-2
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrlab1

/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yychar1 = YYTRANSLATE (yychar);				\
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { 								\
      yyerror ("syntax error: cannot back up");			\
      YYERROR;							\
    }								\
while (0)

#define YYTERROR	1
#define YYERRCODE	256

/* YYLLOC_DEFAULT -- Compute the default location (before the actions
   are run).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)           \
  Current.first_line   = Rhs[1].first_line;      \
  Current.first_column = Rhs[1].first_column;    \
  Current.last_line    = Rhs[N].last_line;       \
  Current.last_column  = Rhs[N].last_column;
#endif

/* YYLEX -- calling `yylex' with the right arguments.  */

#define YYLEX	yylex ()

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (0)
# define YYDSYMPRINT(Args)			\
do {						\
  if (yydebug)					\
    yysymprint Args;				\
} while (0)
/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YYDSYMPRINT(Args)
#endif /* !YYDEBUG */

/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   SIZE_MAX < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#if YYMAXDEPTH == 0
# undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined (__GLIBC__) && defined (_STRING_H)
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
#   if defined (__STDC__) || defined (__cplusplus)
yystrlen (const char *yystr)
#   else
yystrlen (yystr)
     const char *yystr;
#   endif
{
  register const char *yys = yystr;

  while (*yys++ != '\0')
    continue;

  return yys - yystr - 1;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined (__GLIBC__) && defined (_STRING_H) && defined (_GNU_SOURCE)
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
#   if defined (__STDC__) || defined (__cplusplus)
yystpcpy (char *yydest, const char *yysrc)
#   else
yystpcpy (yydest, yysrc)
     char *yydest;
     const char *yysrc;
#   endif
{
  register char *yyd = yydest;
  register const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

#endif /* !YYERROR_VERBOSE */



#if YYDEBUG
/*-----------------------------.
| Print this symbol on YYOUT.  |
`-----------------------------*/

static void
#if defined (__STDC__) || defined (__cplusplus)
yysymprint (FILE* yyout, int yytype, YYSTYPE yyvalue)
#else
yysymprint (yyout, yytype, yyvalue)
    FILE* yyout;
    int yytype;
    YYSTYPE yyvalue;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvalue;

  if (yytype < YYNTOKENS)
    {
      YYFPRINTF (yyout, "token %s (", yytname[yytype]);
# ifdef YYPRINT
      YYPRINT (yyout, yytoknum[yytype], yyvalue);
# endif
    }
  else
    YYFPRINTF (yyout, "nterm %s (", yytname[yytype]);

  switch (yytype)
    {
      default:
        break;
    }
  YYFPRINTF (yyout, ")");
}
#endif /* YYDEBUG. */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
#if defined (__STDC__) || defined (__cplusplus)
yydestruct (int yytype, YYSTYPE yyvalue)
#else
yydestruct (yytype, yyvalue)
    int yytype;
    YYSTYPE yyvalue;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvalue;

  switch (yytype)
    {
      default:
        break;
    }
}



/* The user can define YYPARSE_PARAM as the name of an argument to be passed
   into yyparse.  The argument should have type void *.
   It should actually point to an object.
   Grammar actions can access the variable by casting it
   to the proper pointer type.  */

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
#  define YYPARSE_PARAM_ARG void *YYPARSE_PARAM
#  define YYPARSE_PARAM_DECL
# else
#  define YYPARSE_PARAM_ARG YYPARSE_PARAM
#  define YYPARSE_PARAM_DECL void *YYPARSE_PARAM;
# endif
#else /* !YYPARSE_PARAM */
# define YYPARSE_PARAM_ARG
# define YYPARSE_PARAM_DECL
#endif /* !YYPARSE_PARAM */

/* Prevent warning if -Wstrict-prototypes.  */
#ifdef __GNUC__
# ifdef YYPARSE_PARAM
int yyparse (void *);
# else
int yyparse (void);
# endif
#endif


/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of parse errors so far.  */
int yynerrs;


int
yyparse (YYPARSE_PARAM_ARG)
     YYPARSE_PARAM_DECL
{
  
  register int yystate;
  register int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Lookahead token as an internal (translated) token number.  */
  int yychar1 = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  short	yyssa[YYINITDEPTH];
  short *yyss = yyssa;
  register short *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  register YYSTYPE *yyvsp;



#define YYPOPSTACK   (yyvsp--, yyssp--)

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* When reducing, the number of symbols on the RHS of the reduced
     rule.  */
  int yylen;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed. so pushing a state here evens the stacks.
     */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyssp >= yyss + yystacksize - 1)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack. Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	short *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyoverflowlab;
# else
      /* Extend the stack our own way.  */
      if (yystacksize >= YYMAXDEPTH)
	goto yyoverflowlab;
      yystacksize *= 2;
      if (yystacksize > YYMAXDEPTH)
	yystacksize = YYMAXDEPTH;

      {
	short *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyoverflowlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyssp >= yyss + yystacksize - 1)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* yychar is either YYEMPTY or YYEOF
     or a valid token in external form.  */

  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  /* Convert token to internal form (in yychar1) for indexing tables with.  */

  if (yychar <= 0)		/* This means end of input.  */
    {
      yychar1 = 0;
      yychar = YYEOF;		/* Don't call YYLEX any more.  */

      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yychar1 = YYTRANSLATE (yychar);

      /* We have to keep this `#if YYDEBUG', since we use variables
	 which are defined only if `YYDEBUG' is set.  */
      YYDPRINTF ((stderr, "Next token is "));
      YYDSYMPRINT ((stderr, yychar1, yylval));
      YYDPRINTF ((stderr, "\n"));
    }

  /* If the proper action on seeing token YYCHAR1 is to reduce or to
     detect an error, take that action.  */
  yyn += yychar1;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yychar1)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */
  YYDPRINTF ((stderr, "Shifting token %d (%s), ",
	      yychar, yytname[yychar1]));

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;


  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  yystate = yyn;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];



#if YYDEBUG
  /* We have to keep this `#if YYDEBUG', since we use variables which
     are defined only if `YYDEBUG' is set.  */
  if (yydebug)
    {
      int yyi;

      YYFPRINTF (stderr, "Reducing via rule %d (line %d), ",
		 yyn - 1, yyrline[yyn]);

      /* Print the symbols being reduced, and their result.  */
      for (yyi = yyprhs[yyn]; yyrhs[yyi] >= 0; yyi++)
	YYFPRINTF (stderr, "%s ", yytname[yyrhs[yyi]]);
      YYFPRINTF (stderr, " -> %s\n", yytname[yyr1[yyn]]);
    }
#endif
  switch (yyn)
    {
        case 2:
    {
                Vsh_trace("Yacc","inputunit::1");
                global_command = yyvsp[-1].command;
                if (global_command != VNULL) {
                    if (global_command->type == cm_simple) {
                        wordTmp = yyvsp[-1].command->value.Simple->words;
                        memset(buf, VNULL_SYMBOL, sizeof(buf));
                        while (wordTmp != VNULL) {
                            wTmp = wordTmp->word;
                            if (wTmp != VNULL) {
                                strcat(buf,wTmp->word); 
                                wordTmp = wordTmp->next;
                            } else {
                                wordTmp = VNULL;
                            }
                            if (wordTmp != VNULL) strcat(buf," "); 
                        }
                        Vnm_print(1,"Vsh: will execv: <%s>\n", buf);
                        Vsh_parseHandoff(buf);
                        Vsh_yyexecute(global_command);
                        dispose_command (global_command);
                    }
                }
                YYACCEPT;
            }
    break;

  case 3:
    {
                Vsh_trace("Yacc","inputunit::2");
                global_command = (COMMAND*)VNULL;
                YYACCEPT;
            }
    break;

  case 4:
    {
                Vsh_trace("Yacc","inputunit::3");
                global_command = (COMMAND*)VNULL;
                cmdKey = 2;
                YYACCEPT;
            }
    break;

  case 5:
    {
                Vsh_trace("Yacc","inputunit::4");
                global_command = (COMMAND*)VNULL;
                YYABORT;
            }
    break;

  case 6:
    {
                  Vsh_trace("Yacc","simple_list::1");
                  yyval.command = yyvsp[0].command;
              }
    break;

  case 7:
    {
                  Vsh_trace("Yacc","simple_list::2");
                  yyval.command = yyvsp[-1].command;
              }
    break;

  case 8:
    {
                  Vsh_trace("Yacc","simple_list::3");
                  yyval.command = yyvsp[-1].command;
              }
    break;

  case 9:
    {
                   Vsh_trace("Yacc","simple_list1::1");
                   yyval.command = yyvsp[-3].command;
               }
    break;

  case 10:
    {
                   Vsh_trace("Yacc","simple_list1::2");
                   yyval.command = yyvsp[-3].command;
               }
    break;

  case 11:
    {
                   Vsh_trace("Yacc","simple_list1::3");
                   yyval.command = yyvsp[-2].command;
               }
    break;

  case 12:
    {
                   Vsh_trace("Yacc","simple_list1::4");
                   yyval.command = yyvsp[-2].command;
               }
    break;

  case 13:
    {
                   Vsh_trace("Yacc","simple_list1::5");
                   yyval.command = yyvsp[0].command;
               }
    break;

  case 14:
    {
                   Vsh_trace("Yacc","simple_list1::6");
                   yyval.command = yyvsp[0].command;
               }
    break;

  case 15:
    {
               Vsh_trace("Yacc","pipeline::1");
               yyval.command = yyvsp[-3].command;
           }
    break;

  case 16:
    {
               Vsh_trace("Yacc","pipeline::2");
               yyval.command = yyvsp[0].command;
           }
    break;

  case 17:
    {
           Vsh_trace("Yacc","list::1");
           yyval.command = yyvsp[0].command;
       }
    break;

  case 18:
    {
            Vsh_trace("Yacc","list0::1");
            yyval.command = yyvsp[0].command;
        }
    break;

  case 19:
    {
            Vsh_trace("Yacc","list0::2");
            yyval.command = yyvsp[-2].command;
        }
    break;

  case 20:
    {
            Vsh_trace("Yacc","list0::3");
            yyval.command = yyvsp[-2].command;
        }
    break;

  case 21:
    {
            Vsh_trace("Yacc","list0::4");
            yyval.command = yyvsp[-2].command;
        }
    break;

  case 22:
    {
            Vsh_trace("Yacc","list1::1");
            yyval.command = yyvsp[-3].command;
        }
    break;

  case 23:
    {
            Vsh_trace("Yacc","list1::2");
            yyval.command = yyvsp[-3].command;
        }
    break;

  case 24:
    {
            Vsh_trace("Yacc","list1::3");
            yyval.command = yyvsp[-3].command;
        }
    break;

  case 25:
    {
            Vsh_trace("Yacc","list1::4");
            yyval.command = yyvsp[-3].command;
        }
    break;

  case 26:
    {
            Vsh_trace("Yacc","list1::5");
            yyval.command = yyvsp[-3].command;
        }
    break;

  case 27:
    {
            Vsh_trace("Yacc","list1::6");
            yyval.command = yyvsp[0].command;
        }
    break;

  case 28:
    {
            Vsh_trace("Yacc","list1::7");
            yyval.command = yyvsp[0].command;
        }
    break;

  case 29:
    {
                  Vsh_trace("Yacc","redirection::1");
                  redir.filename = yyvsp[0].word;
                  yyval.redirect = make_redirection(1, r_output_direction, redir);
              }
    break;

  case 30:
    {
                  Vsh_trace("Yacc","redirection::2");
                  redir.filename = yyvsp[0].word;
                  yyval.redirect = make_redirection(1, r_appending_to, redir);
              }
    break;

  case 31:
    {
                   Vsh_trace("Yacc","redirections::1");
                   yyval.redirect = yyvsp[0].redirect;
               }
    break;

  case 32:
    {
                   Vsh_trace("Yacc","redirections::2");
                   yyval.redirect = yyvsp[-1].redirect;
               }
    break;

  case 33:
    {
              Vsh_trace("Yacc","command::1");
              yyval.command = clean_simple_command(yyvsp[0].command);
          }
    break;

  case 34:
    {
              Vsh_trace("Yacc","command::2");
              yyval.command = yyvsp[0].command;
          }
    break;

  case 35:
    {
                     Vsh_trace("Yacc","simple_command::1");
                     yyval.command = make_simple_command(yyvsp[0].element,(COMMAND*)VNULL);
                 }
    break;

  case 36:
    {
                     Vsh_trace("Yacc","simple_command::2");
                     yyval.command = make_simple_command(yyvsp[0].element,yyvsp[-1].command);
                 }
    break;

  case 37:
    {
                    Vsh_trace("Yacc","shell_command::1");
                    yyval.command = yyvsp[0].command;
                }
    break;

  case 38:
    {
                    Vsh_trace("Yacc","shell_command::2");
                    yyval.command = yyvsp[-1].command;
                }
    break;

  case 39:
    {
                      Vsh_trace("Yacc","shell_command_1::1");
                      yyval.command = VNULL;
                  }
    break;

  case 40:
    {
                      Vsh_trace("Yacc","shell_command_1::2");
                      yyval.command = VNULL;
                  }
    break;

  case 41:
    {
                      Vsh_trace("Yacc","shell_command_1::3");
                      yyval.command = VNULL;
                  }
    break;

  case 42:
    {
                      Vsh_trace("Yacc","shell_command_1::4");
                      yyval.command = VNULL;
                  }
    break;

  case 43:
    {
                      Vsh_trace("Yacc","shell_command_1::5");
                      yyval.command = VNULL;
                  }
    break;

  case 44:
    {
                      Vsh_trace("Yacc","shell_command_1::6");
                      yyval.command = VNULL;
                  }
    break;

  case 45:
    {
                      Vsh_trace("Yacc","shell_command_1::7");
                      yyval.command = VNULL;
                  }
    break;

  case 46:
    {
                      Vsh_trace("Yacc","shell_command_1::8");
                      yyval.command = VNULL;
                  }
    break;

  case 47:
    {
                             Vsh_trace("Yacc","simple_command_element::1");
                             yyval.element.word     = yyvsp[0].word;
                             yyval.element.redirect = 0;
                         }
    break;

  case 48:
    {
                             Vsh_trace("Yacc","simple_command_element::2");
                             yyval.element.word     = yyvsp[0].word;
                             yyval.element.redirect = 0;
                         }
    break;

  case 49:
    {
                             Vsh_trace("Yacc","simple_command_element::3");
                             yyval.element.word     = 0;
                             yyval.element.redirect = yyvsp[0].redirect;
                         }
    break;

  case 50:
    {
               Vsh_trace("Yacc","subshell::1");
               yyval.command = yyvsp[-1].command;
           }
    break;

  case 51:
    {
                    Vsh_trace("Yacc","group_command::1");
                    yyval.command = yyvsp[-1].command;
                }
    break;

  case 52:
    {
                   Vsh_trace("Yacc","function_def::1");
                   yyval.command = VNULL;
               }
    break;

  case 53:
    {
                   Vsh_trace("Yacc","function_def::2");
                   yyval.command = VNULL;
               }
    break;

  case 54:
    {
                 Vsh_trace("Yacc","if_command::1");
                 yyval.command = VNULL;
             }
    break;

  case 55:
    {
                 Vsh_trace("Yacc","if_command::2");
                 yyval.command = VNULL;
             }
    break;

  case 56:
    {
                 Vsh_trace("Yacc","if_command::3");
                 yyval.command = VNULL;
             }
    break;

  case 57:
    {
                  Vsh_trace("Yacc","elif_clause::1");
                  yyval.command = VNULL;
              }
    break;

  case 58:
    {
                  Vsh_trace("Yacc","elif_clause::2");
                  yyval.command = VNULL;
              }
    break;

  case 59:
    {
                  Vsh_trace("Yacc","elif_clause::3");
                  yyval.command = VNULL;
              }
    break;

  case 60:
    {
                   Vsh_trace("Yacc","case_command::1");
                   yyval.command = VNULL;
               }
    break;

  case 61:
    {
                   Vsh_trace("Yacc","case_command::2");
                   yyval.command = VNULL;
               }
    break;

  case 62:
    {
                   Vsh_trace("Yacc","case_command::3");
                   yyval.command = VNULL;
               }
    break;

  case 63:
    {
                    Vsh_trace("Yacc","case_clause_1::1");
                    yyval.pattern = VNULL;
                }
    break;

  case 64:
    {
                    Vsh_trace("Yacc","case_clause_1::2");
                    yyval.pattern = VNULL;
                }
    break;

  case 65:
    {
                     Vsh_trace("Yacc","pattern_list_1::1");
                     yyval.pattern = VNULL;
                 }
    break;

  case 66:
    {
                     Vsh_trace("Yacc","pattern_list_1::2");
                     yyval.pattern = VNULL;
                 }
    break;

  case 67:
    {
                     Vsh_trace("Yacc","pattern_list_1::3");
                     yyval.pattern = VNULL;
                 }
    break;

  case 68:
    {
                     Vsh_trace("Yacc","pattern_list_1::4");
                     yyval.pattern = VNULL;
                 }
    break;

  case 69:
    {
                           Vsh_trace("Yacc","case_clause_sequence::1");
                           yyval.pattern = VNULL;
                       }
    break;

  case 70:
    {
                           Vsh_trace("Yacc","case_clause_sequence::2");
                           yyval.pattern = VNULL;
                       }
    break;

  case 71:
    {
                   Vsh_trace("Yacc","pattern_list::1");
                   yyval.pattern = VNULL;
               }
    break;

  case 72:
    {
                   Vsh_trace("Yacc","pattern_list::2");
                   yyval.pattern = VNULL;
               }
    break;

  case 73:
    {
                   Vsh_trace("Yacc","pattern_list::3");
                   yyval.pattern = VNULL;
               }
    break;

  case 74:
    {
                   Vsh_trace("Yacc","pattern_list::4");
                   yyval.pattern = VNULL;
               }
    break;

  case 75:
    {
              Vsh_trace("Yacc","pattern::1");
              yyval.word_list = VNULL;
          }
    break;

  case 76:
    {
              Vsh_trace("Yacc","pattern::2");
              yyval.word_list = VNULL;
          }
    break;

  case 77:
    {
                  Vsh_trace("Yacc","for_command::1");
                  yyval.command = VNULL;
              }
    break;

  case 78:
    {
                  Vsh_trace("Yacc","for_command::2");
                  yyval.command = VNULL;
              }
    break;

  case 79:
    {
                  Vsh_trace("Yacc","for_command::3");
                  yyval.command = VNULL;
              }
    break;

  case 80:
    {
                  Vsh_trace("Yacc","for_command::4");
                  yyval.command = VNULL;
              }
    break;

  case 81:
    {
                  Vsh_trace("Yacc","for_command::5");
                  yyval.command = VNULL;
              }
    break;

  case 82:
    {
                  Vsh_trace("Yacc","for_command::6");
                  yyval.command = VNULL;
              }
    break;

  case 83:
    {
                    Vsh_trace("Yacc","while_command::1");
                    yyval.command = VNULL;
                }
    break;

  case 84:
    {
                    Vsh_trace("Yacc","until_command::1");
                    yyval.command = VNULL;
                }
    break;

  case 85:
    {
            Vsh_trace("Yacc","words::1");
            yyval.word_list = (WORD_LIST*)VNULL;
        }
    break;

  case 86:
    {
            Vsh_trace("Yacc","words::2");
            yyval.word_list = make_word_list(yyvsp[0].word,yyvsp[-1].word_list);
        }
    break;

  case 87:
    {
           Vsh_trace("Yacc","term::1");
       }
    break;

  case 88:
    {
           Vsh_trace("Yacc","term::2");
       }
    break;

  case 89:
    {
           Vsh_trace("Yacc","term::3");
       }
    break;

  case 90:
    {
              Vsh_trace("Yacc","newlines::1");
          }
    break;

  case 91:
    {
              Vsh_trace("Yacc","newlines::2");
          }
    break;


    }

/* Line 1016 of /usr/share/bison/yacc.c.  */
#line 2236 "y.tab.c"

  yyvsp -= yylen;
  yyssp -= yylen;


#if YYDEBUG
  if (yydebug)
    {
      short *yyssp1 = yyss - 1;
      YYFPRINTF (stderr, "state stack now");
      while (yyssp1 != yyssp)
	YYFPRINTF (stderr, " %d", *++yyssp1);
      YYFPRINTF (stderr, "\n");
    }
#endif

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (YYPACT_NINF < yyn && yyn < YYLAST)
	{
	  YYSIZE_T yysize = 0;
	  int yytype = YYTRANSLATE (yychar);
	  char *yymsg;
	  int yyx, yycount;

	  yycount = 0;
	  /* Start YYX at -YYN if negative to avoid negative indexes in
	     YYCHECK.  */
	  for (yyx = yyn < 0 ? -yyn : 0;
	       yyx < (int) (sizeof (yytname) / sizeof (char *)); yyx++)
	    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	      yysize += yystrlen (yytname[yyx]) + 15, yycount++;
	  yysize += yystrlen ("parse error, unexpected ") + 1;
	  yysize += yystrlen (yytname[yytype]);
	  yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg != 0)
	    {
	      char *yyp = yystpcpy (yymsg, "parse error, unexpected ");
	      yyp = yystpcpy (yyp, yytname[yytype]);

	      if (yycount < 5)
		{
		  yycount = 0;
		  for (yyx = yyn < 0 ? -yyn : 0;
		       yyx < (int) (sizeof (yytname) / sizeof (char *));
		       yyx++)
		    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
		      {
			const char *yyq = ! yycount ? ", expecting " : " or ";
			yyp = yystpcpy (yyp, yyq);
			yyp = yystpcpy (yyp, yytname[yyx]);
			yycount++;
		      }
		}
	      yyerror (yymsg);
	      YYSTACK_FREE (yymsg);
	    }
	  else
	    yyerror ("parse error; also virtual memory exhausted");
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror ("parse error");
    }
  goto yyerrlab1;


/*----------------------------------------------------.
| yyerrlab1 -- error raised explicitly by an action.  |
`----------------------------------------------------*/
yyerrlab1:
  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      /* Return failure if at end of input.  */
      if (yychar == YYEOF)
        {
	  /* Pop the error token.  */
          YYPOPSTACK;
	  /* Pop the rest of the stack.  */
	  while (yyssp > yyss)
	    {
	      YYDPRINTF ((stderr, "Error: popping "));
	      YYDSYMPRINT ((stderr,
			    yystos[*yyssp],
			    *yyvsp));
	      YYDPRINTF ((stderr, "\n"));
	      yydestruct (yystos[*yyssp], *yyvsp);
	      YYPOPSTACK;
	    }
	  YYABORT;
        }

      YYDPRINTF ((stderr, "Discarding token %d (%s).\n",
		  yychar, yytname[yychar1]));
      yydestruct (yychar1, yylval);
      yychar = YYEMPTY;
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */

  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;

      YYDPRINTF ((stderr, "Error: popping "));
      YYDSYMPRINT ((stderr,
		    yystos[*yyssp], *yyvsp));
      YYDPRINTF ((stderr, "\n"));

      yydestruct (yystos[yystate], *yyvsp);
      yyvsp--;
      yystate = *--yyssp;


#if YYDEBUG
      if (yydebug)
	{
	  short *yyssp1 = yyss - 1;
	  YYFPRINTF (stderr, "Error: state stack now");
	  while (yyssp1 != yyssp)
	    YYFPRINTF (stderr, " %d", *++yyssp1);
	  YYFPRINTF (stderr, "\n");
	}
#endif
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  YYDPRINTF ((stderr, "Shifting error token, "));

  *++yyvsp = yylval;


  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*----------------------------------------------.
| yyoverflowlab -- parser overflow comes here.  |
`----------------------------------------------*/
yyoverflowlab:
  yyerror ("parser stack overflow");
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}





