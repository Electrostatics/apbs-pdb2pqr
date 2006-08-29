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

#ifndef BISON_Y_TAB_H
# define BISON_Y_TAB_H

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
/* Line 1281 of /usr/share/bison/yacc.c.  */
#line 114 "y.tab.h"
# define YYSTYPE yystype
#endif

extern YYSTYPE yylval;


#endif /* not BISON_Y_TAB_H */

