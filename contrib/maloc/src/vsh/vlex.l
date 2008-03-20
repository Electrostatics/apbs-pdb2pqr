/*
 * ***************************************************************************
 * MALOC = < Minimal Abstraction Layer for Object-oriented C >
 * Copyright (C) 1994--2008 Michael Holst
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * rcsid="$Id: vlex.l,v 1.7 2008/03/12 05:13:58 fetk Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     vlex.l
 *
 * Purpose:  Class Vsh: LEX Specification
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

UND [_]
LET [_A-Za-z]
DIG [0-9]
ULD [_A-Za-z0-9]
SYM [-+$*.:/\\{}~^|@<>%#]

%{

#include "vsh_p.h"

VEMBED(rcsid="$Id: vlex.l,v 1.7 2008/03/12 05:13:58 fetk Exp $")

#define YY_NO_UNPUT 1
#undef YY_INPUT
#define YY_INPUT(buf,result,max_size) VSH_INPUT(buf,result,max_size)

#include "vyacc.h"

%}

/*
 * ***************************************************************************
 * Initial states
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * LEX Analyzer rules
 * ***************************************************************************
 */

%%

    /*
     * ***********************************************************************
     * Whitespace
     * ***********************************************************************
     */

[ \t]+ { /* no-op: ignore white space between tokens */ }

^[ ]*\n { Vsh_trace("Lex","line of zero or more blanks with a newline"); }

\\\n { Vsh_trace("Lex","an escaped newline"); }

    /*
     * ***********************************************************************
     * Comments
     * ***********************************************************************
     */

    /* ignore shell-style "to-end-of-line" comments (OReilly, page 31) */
#.*\n { Vsh_trace("Lex","single-line ``#'' comment"); }

    /*
     * ***********************************************************************
     * Legal 2-character symbol tokens (must encode them...)
     * ***********************************************************************
     */

";;" {
    Vsh_trace("Lex",yytext);
    return SEMI_SEMI;
}
"&&" {
    Vsh_trace("Lex",yytext);
    return AND_AND;
}
"||" {
    Vsh_trace("Lex",yytext);
    return OR_OR;
}
"<<" {
    Vsh_trace("Lex",yytext);
    return LESS_LESS;
}
">>" {
    Vsh_trace("Lex",yytext);
    return GREATER_GREATER;
}
"<&" {
    Vsh_trace("Lex",yytext);
    return LESS_AND;
}
">&" {
    Vsh_trace("Lex",yytext);
    return GREATER_AND;
}
"&>" {
    Vsh_trace("Lex",yytext);
    return AND_GREATER;
}
"<>" {
    Vsh_trace("Lex",yytext);
    return LESS_GREATER;
}
">|" {
    Vsh_trace("Lex",yytext);
    return GREATER_BAR;
}

    /*
     * ***********************************************************************
     * Legal 3-character symbol tokens (must encode them...)
     * ***********************************************************************
     */

"<<-" {
    Vsh_trace("Lex",yytext);
    return LESS_LESS_MINUS;
}

    /*
     * ***********************************************************************
     * Legal 1-character symbol tokens (just forward them to the parser...)
     * ***********************************************************************
     */

"="  |
"("  |
")"  |
"<"  |
">"  |
"|"  |
"&"  |
";"  |
"\n" {
    Vsh_trace("Lex",yytext);
    return yytext[0];
}

    /*
     * ***********************************************************************
     * Reserved words
     * ***********************************************************************
     */

"if" {
    Vsh_trace("Lex","IF");
    return IF;
}
"then" {
    Vsh_trace("Lex","THEN");
    return THEN;
}
"else" {
    Vsh_trace("Lex","ELSE");
    return ELSE;
}
"elif" {
    Vsh_trace("Lex","ELIF");
    return ELIF;
}
"fi" {
    Vsh_trace("Lex","FI");
    return FI;
}
"case" {
    Vsh_trace("Lex","CASE");
    return CASE;
}
"in" {
    Vsh_trace("Lex","IN");
    return IN;
}
"esac" {
    Vsh_trace("Lex","ESAC");
    return ESAC;
}
"for" {
    Vsh_trace("Lex","FOR");
    return FOR;
}
"while" {
    Vsh_trace("Lex","WHILE");
    return WHILE;
}
"until" {
    Vsh_trace("Lex","UNTIL");
    return UNTIL;
}
"do" {
    Vsh_trace("Lex","DO");
    return DO;
}
"done" {
    Vsh_trace("Lex","DONE");
    return DONE;
}
"{" {
    Vsh_trace("Lex",yytext);
    return yytext[0];
}
"}" {
    Vsh_trace("Lex",yytext);
    return yytext[0];
}

    /*
     * ***********************************************************************
     * Tokens
     * ***********************************************************************
     */

\'([^']|\\\')*\' {
    Vsh_trace("Lex","WORD-SSTRING");
    yylval.word = (WORD_DESC*)malloc(sizeof(WORD_DESC));
    yylval.word->word = (char*)malloc(1024*sizeof(char));
    strcpy(yylval.word->word,yytext);
    yylval.word->dollar_present = 0;
    yylval.word->quoted         = 1;
    yylval.word->assignment     = 0;
    return WORD;
}

\"([^"]|\\\")*\" {
    Vsh_trace("Lex","WORD-DSTRING");
    yylval.word = (WORD_DESC*)malloc(sizeof(WORD_DESC));
    yylval.word->word = (char*)malloc(1024*sizeof(char));
    strcpy(yylval.word->word,yytext);
    yylval.word->dollar_present = 0;
    yylval.word->quoted         = 1;
    yylval.word->assignment     = 0;
    return WORD;
}

\`([^`]|\\\`)*\` {
    Vsh_trace("Lex","WORD-EXEC");
    yylval.word = (WORD_DESC*)malloc(sizeof(WORD_DESC));
    yylval.word->word = (char*)malloc(1024*sizeof(char));
    strcpy(yylval.word->word,yytext);
    yylval.word->dollar_present = 0;
    yylval.word->quoted         = 1;
    yylval.word->assignment     = 0;
    return WORD;
}

"["([^[]|\\\[)*"]" {
    Vsh_trace("Lex","WORD-CONDITION");
    yylval.word = (WORD_DESC*)malloc(sizeof(WORD_DESC));
    yylval.word->word = (char*)malloc(1024*sizeof(char));
    strcpy(yylval.word->word,yytext);
    yylval.word->dollar_present = 0;
    yylval.word->quoted         = 0;
    yylval.word->assignment     = 0;
    return WORD;
}

{LET}+{ULD}*[=] {
    /* struct symtab *sp = symlook(yytext); */
    /* yylval.symp = sp; */
    Vsh_trace("Lex","WORD-ASSIGN");
    yylval.word = (WORD_DESC*)malloc(sizeof(WORD_DESC));
    yylval.word->word = (char*)malloc(1024*sizeof(char));
    strcpy(yylval.word->word,yytext);
    VASSERT( yylval.word->word[strlen(yylval.word->word)-1] == '=' );
    yylval.word->word[strlen(yylval.word->word)-1] = '\0';
    yylval.word->dollar_present = 0;
    yylval.word->quoted         = 0;
    yylval.word->assignment     = 1;
    return ASSIGNMENT_WORD;
}

({SYM}|{ULD})+ {
    /* struct symtab *sp = symlook(yytext); */
    /* yylval.symp = sp; */
    Vsh_trace("Lex","WORD");
    yylval.word = (WORD_DESC*)malloc(sizeof(WORD_DESC));
    yylval.word->word = (char*)malloc(1024*sizeof(char));
    strcpy(yylval.word->word,yytext);
    yylval.word->dollar_present = 0;
    yylval.word->quoted         = 0;
    yylval.word->assignment     = 0;
    return WORD;
}

<<EOF>> {
    Vsh_trace("Lex","vsh_EOF");
    return vsh_EOF;
}

    /*
     * ***********************************************************************
     * Catchall (error exit)
     * ***********************************************************************
     */

. {
    yyerror("bad token");
}

%%

