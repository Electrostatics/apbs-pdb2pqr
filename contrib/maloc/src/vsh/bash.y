/* Yacc grammar for bash. */

/* Copyright (C) 1989 Free Software Foundation, Inc.

   This file is part of GNU Bash, the Bourne Again SHell.

   Bash is free software; you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free
   Software Foundation; either version 1, or (at your option) any later
   version.

   Bash is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   for more details.

   You should have received a copy of the GNU General Public License along
   with Bash; see the file LICENSE.  If not, write to the Free Software
   Foundation, 675 Mass Ave, Cambridge, MA 02139, USA. */

%{
#include <stdio.h>
#include "bashtypes.h"
#include <signal.h>
#include "bashansi.h"
#include "shell.h"
#include "flags.h"
#include "input.h"

#if defined (READLINE)
#  include <readline/readline.h>
#endif /* READLINE */

#if defined (HISTORY)
#  include "bashhist.h"
#  include <readline/history.h>
#endif /* HISTORY */

#if defined (JOB_CONTROL)
#  include "jobs.h"
#endif /* JOB_CONTROL */

#if defined (ALIAS)
#  include "alias.h"
#endif /* ALIAS */

#if defined (PROMPT_STRING_DECODE)
#include <sys/param.h>
#include <time.h>
#include "maxpath.h"
#endif /* PROMPT_STRING_DECODE */

#define YYDEBUG 1
extern int eof_encountered;
extern int no_line_editing;
extern int current_command_number;
extern int interactive, interactive_shell, login_shell;
extern int posixly_correct;
extern int last_command_exit_value;
extern int interrupt_immediately;
extern char *shell_name, *current_host_name;
extern Function *last_shell_builtin, *this_shell_builtin;
#if defined (READLINE)
extern int bash_readline_initialized;
#endif
#if defined (BUFFERED_INPUT)
extern int bash_input_fd_changed;
#endif

/* **************************************************************** */
/*								    */
/*		    "Forward" declarations			    */
/*								    */
/* **************************************************************** */

/* This is kind of sickening.  In order to let these variables be seen by
   all the functions that need them, I am forced to place their declarations
   far away from the place where they should logically be found. */

static int reserved_word_acceptable ();
static int read_token ();

static void report_syntax_error ();
static void handle_eof_input_unit ();
static void prompt_again ();
static void reset_readline_prompt ();
static void print_prompt ();

/* PROMPT_STRING_POINTER points to one of these, never to an actual string. */
char *ps1_prompt, *ps2_prompt;

/* Handle on the current prompt string.  Indirectly points through
   ps1_ or ps2_prompt. */
char **prompt_string_pointer = (char **)NULL;
char *current_prompt_string;

/* The decoded prompt string.  Used if READLINE is not defined or if
   editing is turned off.  Analogous to current_readline_prompt. */
static char *current_decoded_prompt;

/* The number of lines read from input while creating the current command. */
int current_command_line_count = 0;

/* Variables to manage the task of reading here documents, because we need to
   defer the reading until after a complete command has been collected. */
static REDIRECT *redir_stack[10];
int need_here_doc = 0;

/* Where shell input comes from.  History expansion is performed on each
   line when the shell is interactive. */
static char *shell_input_line = (char *)NULL;
static int shell_input_line_index = 0;
static int shell_input_line_size = 0;	/* Amount allocated for shell_input_line. */
static int shell_input_line_len = 0;	/* strlen (shell_input_line) */

/* Either zero or EOF. */
static int shell_input_line_terminator = 0;

static REDIRECTEE redir;
%}

%union {
  WORD_DESC *word;		/* the word that we read. */
  int number;			/* the number that we read. */
  WORD_LIST *word_list;
  COMMAND *command;
  REDIRECT *redirect;
  ELEMENT element;
  PATTERN_LIST *pattern;
}

/* Reserved words.  Members of the first group are only recognized
   in the case that they are preceded by a list_terminator.  Members
   of the second group are recognized only under special circumstances. */
%token IF THEN ELSE ELIF FI CASE ESAC FOR SELECT WHILE UNTIL DO DONE FUNCTION
%token IN BANG

/* More general tokens. yylex () knows how to make these. */
%token <word> WORD ASSIGNMENT_WORD
%token <number> NUMBER
%token AND_AND OR_OR GREATER_GREATER LESS_LESS LESS_AND
%token GREATER_AND SEMI_SEMI LESS_LESS_MINUS AND_GREATER LESS_GREATER
%token GREATER_BAR

/* The types that the various syntactical units return. */

%type <command> inputunit command pipeline
%type <command> list list0 list1 simple_list simple_list1
%type <command> simple_command shell_command_1 shell_command select_command
%type <command> group_command function_def if_command elif_clause subshell
%type <redirect> redirection redirections
%type <element> simple_command_element
%type <word_list> words pattern 
%type <pattern> pattern_list case_clause_sequence case_clause_1 pattern_list_1

%start inputunit

%left '&' ';' '\n' yacc_EOF
%left AND_AND OR_OR
%right '|'
%%

inputunit:	simple_list '\n'
			{
			  /* Case of regular command.  Discard the error
			     safety net,and return the command just parsed. */
			  global_command = $1;
			  eof_encountered = 0;
			  discard_parser_constructs (0);
			  YYACCEPT;
			}
	|	'\n'
			{
			  /* Case of regular command, but not a very
			     interesting one.  Return a NULL command. */
			  global_command = (COMMAND *)NULL;
			  YYACCEPT;
			}
	|
		error '\n'
			{
			  /* Error during parsing.  Return NULL command. */
			  global_command = (COMMAND *)NULL;
			  eof_encountered = 0;
			  discard_parser_constructs (1);
			  if (interactive)
			    {
			      YYACCEPT;
			    }
			  else
			    {
			      YYABORT;
			    }
			}
	|	yacc_EOF
			{
			  /* Case of EOF seen by itself.  Do ignoreeof or 
			     not. */
			  global_command = (COMMAND *)NULL;
			  handle_eof_input_unit ();
			  YYACCEPT;
			}
	;

words:	
			{ $$ = (WORD_LIST *)NULL; }
	|	words WORD
			{ $$ = make_word_list ($2, $1); }
	;

redirection:	'>' WORD
			{
			  redir.filename = $2;
			  $$ = make_redirection (1, r_output_direction, redir);
			}
	|	'<' WORD
			{
			  redir.filename = $2;
			  $$ = make_redirection (0, r_input_direction, redir);
			}
	|	NUMBER '>' WORD
			{
			  redir.filename = $3;
			  $$ = make_redirection ($1, r_output_direction, redir);
			}
	|	NUMBER '<' WORD
			{
			  redir.filename = $3;
			  $$ = make_redirection ($1, r_input_direction, redir);
			}
	|	GREATER_GREATER WORD
			{
			  redir.filename = $2;
			  $$ = make_redirection (1, r_appending_to, redir);
			}
	|	NUMBER GREATER_GREATER WORD
			{
			  redir.filename = $3;
			  $$ = make_redirection ($1, r_appending_to, redir);
			}
	|	LESS_LESS WORD
			{
			  redir.filename = $2;
			  $$ = make_redirection (0, r_reading_until, redir);
			  redir_stack[need_here_doc++] = $$;
			}
	|	NUMBER LESS_LESS WORD
			{
			  redir.filename = $3;
			  $$ = make_redirection ($1, r_reading_until, redir);
			  redir_stack[need_here_doc++] = $$;
			}
	|	LESS_AND NUMBER
			{
			  redir.dest = $2;
			  $$ = make_redirection (0, r_duplicating_input, redir);
			}
	|	NUMBER LESS_AND NUMBER
			{
			  redir.dest = $3;
			  $$ = make_redirection ($1, r_duplicating_input, redir);
			}
	|	GREATER_AND NUMBER
			{
			  redir.dest = $2;
			  $$ = make_redirection (1, r_duplicating_output, redir);
			}
	|	NUMBER GREATER_AND NUMBER
			{
			  redir.dest = $3;
			  $$ = make_redirection ($1, r_duplicating_output, redir);
			}
	|	LESS_AND WORD
			{
			  redir.filename = $2;
			  $$ = make_redirection (0, r_duplicating_input_word, redir);
			}
	|	NUMBER LESS_AND WORD
			{
			  redir.filename = $3;
			  $$ = make_redirection ($1, r_duplicating_input_word, redir);
			}
	|	GREATER_AND WORD
			{
			  redir.filename = $2;
			  $$ = make_redirection (1, r_duplicating_output_word, redir);
			}
	|	NUMBER GREATER_AND WORD
			{
			  redir.filename = $3;
			  $$ = make_redirection ($1, r_duplicating_output_word, redir);
			}
	|	LESS_LESS_MINUS WORD
			{
			  redir.filename = $2;
			  $$ = make_redirection
			    (0, r_deblank_reading_until, redir);
			  redir_stack[need_here_doc++] = $$;
			}
	|	NUMBER LESS_LESS_MINUS WORD
			{
			  redir.filename = $3;
			  $$ = make_redirection
			    ($1, r_deblank_reading_until, redir);
			  redir_stack[need_here_doc++] = $$;
			}
	|	GREATER_AND '-'
			{
			  redir.dest = 0L;
			  $$ = make_redirection (1, r_close_this, redir);
			}
	|	NUMBER GREATER_AND '-'
			{
			  redir.dest = 0L;
			  $$ = make_redirection ($1, r_close_this, redir);
			}
	|	LESS_AND '-'
			{
			  redir.dest = 0L;
			  $$ = make_redirection (0, r_close_this, redir);
			}
	|	NUMBER LESS_AND '-'
			{
			  redir.dest = 0L;
			  $$ = make_redirection ($1, r_close_this, redir);
			}
	|	AND_GREATER WORD
			{
			  redir.filename = $2;
			  $$ = make_redirection (1, r_err_and_out, redir);
			}
	|	NUMBER LESS_GREATER WORD
			{
			  redir.filename = $3;
			  $$ = make_redirection ($1, r_input_output, redir);
			}
	|	LESS_GREATER WORD
			{
			  REDIRECT *t1, *t2;

			  redir.filename = $2;
			  if (posixly_correct)
			    $$ = make_redirection (0, r_input_output, redir);
			  else
			    {
			      t1 = make_redirection (0, r_input_direction, redir);
			      redir.filename = copy_word ($2);
			      t2 = make_redirection (1, r_output_direction, redir);
			      t1->next = t2;
			      $$ = t1;
			    }
			}			  
	|	GREATER_BAR WORD
			{
			  redir.filename = $2;
			  $$ = make_redirection (1, r_output_force, redir);
			}
	|	NUMBER GREATER_BAR WORD
			{
			  redir.filename = $3;
			  $$ = make_redirection ($1, r_output_force, redir);
			}
	;

simple_command_element: WORD
			{ $$.word = $1; $$.redirect = 0; }
	|	ASSIGNMENT_WORD
			{ $$.word = $1; $$.redirect = 0; }
	|	redirection
			{ $$.redirect = $1; $$.word = 0; }
	;

redirections:	redirection
			{
			  $$ = $1;
			}
	|	redirections redirection
			{ 
			  register REDIRECT *t = $1;

			  while (t->next)
			    t = t->next;
			  t->next = $2; 
			  $$ = $1;
			}
	;

simple_command:	simple_command_element
			{ $$ = make_simple_command ($1, (COMMAND *)NULL); }
	|	simple_command simple_command_element
			{ $$ = make_simple_command ($2, $1); }
	;

command:	simple_command
			{ $$ = clean_simple_command ($1); }
	|	shell_command
			{ $$ = $1; }
	;

shell_command:	shell_command_1
			{ $$ = $1; }
	|	shell_command_1 redirections
			{
			  if ($1->redirects)
			    {
			      register REDIRECT *t;
			      for (t = $1->redirects; t->next; t = t->next)
				;
			      t->next = $2;
			    }
			  else
			    $1->redirects = $2;
			  $$ = $1;
			}
	;

shell_command_1: FOR WORD newlines DO list DONE
			{ $$ = make_for_command ($2, add_string_to_list ("\"$@\"", (WORD_LIST *)NULL), $5); }
	|	FOR WORD newlines '{' list '}'
			{ $$ = make_for_command ($2, add_string_to_list ("$@", (WORD_LIST *)NULL), $5); }
	|	FOR WORD ';' newlines DO list DONE
			{ $$ = make_for_command ($2, add_string_to_list ("\"$@\"", (WORD_LIST *)NULL), $6); }
	|	FOR WORD ';' newlines '{' list '}'
			{ $$ = make_for_command ($2, add_string_to_list ("\"$@\"", (WORD_LIST *)NULL), $6); }
	|	FOR WORD newlines IN words list_terminator newlines DO list DONE
			{ $$ = make_for_command ($2, REVERSE_LIST ($5, WORD_LIST *), $9); }
	|	FOR WORD newlines IN words list_terminator newlines '{' list '}'
			{ $$ = make_for_command ($2, REVERSE_LIST ($5, WORD_LIST *), $9); }

	|	CASE WORD newlines IN newlines ESAC
			{ $$ = make_case_command ($2, (PATTERN_LIST *)NULL); }
	|	CASE WORD newlines IN case_clause_sequence newlines ESAC
			{ $$ = make_case_command ($2, $5); }
	|	CASE WORD newlines IN case_clause_1 ESAC
			{ $$ = make_case_command ($2, $5); }
 	|	WHILE list DO list DONE
			{ $$ = make_while_command ($2, $4); }
	|	UNTIL list DO list DONE
			{ $$ = make_until_command ($2, $4); }
	|	select_command
			{ $$ = $1; }
	|	if_command
			{ $$ = $1; }
	|	subshell
			{ $$ = $1; }
	|	group_command
			{ $$ = $1; }
	|	function_def
			{ $$ = $1; }
	;

select_command:	SELECT WORD newlines DO list DONE
			{
#if defined (SELECT_COMMAND)
			  $$ = make_select_command ($2, add_string_to_list ("\"$@\"", (WORD_LIST *)NULL), $5);
#endif
			}
	|	SELECT WORD newlines '{' list '}'
			{
#if defined (SELECT_COMMAND)
			  $$ = make_select_command ($2, add_string_to_list ("$@", (WORD_LIST *)NULL), $5);
#endif
			}
	|	SELECT WORD ';' newlines DO list DONE
			{
#if defined (SELECT_COMMAND)
			  $$ = make_select_command ($2, add_string_to_list ("\"$@\"", (WORD_LIST *)NULL), $6);
#endif
			}
	|	SELECT WORD ';' newlines '{' list '}'
			{
#if defined (SELECT_COMMAND)
			  $$ = make_select_command ($2, add_string_to_list ("\"$@\"", (WORD_LIST *)NULL), $6);
#endif
			}
	|	SELECT WORD newlines IN words list_terminator newlines DO list DONE
			{
#if defined (SELECT_COMMAND)
			  $$ = make_select_command ($2, (WORD_LIST *)reverse_list ($5), $9);
#endif
			}
	|	SELECT WORD newlines IN words list_terminator newlines '{' list '}'
			{
#if defined (SELECT_COMMAND)
			  $$ = make_select_command ($2, (WORD_LIST *)reverse_list ($5), $9);
#endif
			}
	;

function_def:	WORD '(' ')' newlines group_command
			{ $$ = make_function_def ($1, $5); }

	|	WORD '(' ')' newlines group_command redirections
			{ $5->redirects = $6; $$ = make_function_def ($1, $5); }

	|	FUNCTION WORD '(' ')' newlines group_command
			{ $$ = make_function_def ($2, $6); }

	|	FUNCTION WORD '(' ')' newlines group_command redirections
			{ $6->redirects = $7; $$ = make_function_def ($2, $6); }

	|	FUNCTION WORD newlines group_command
			{ $$ = make_function_def ($2, $4); }

	|	FUNCTION WORD newlines group_command redirections
			{ $4->redirects = $5; $$ = make_function_def ($2, $4); }
	;

subshell:	'(' list ')'
			{ $2->flags |= CMD_WANT_SUBSHELL; $$ = $2; }
	;
	
if_command:	IF list THEN list FI
			{ $$ = make_if_command ($2, $4, (COMMAND *)NULL); }
	|	IF list THEN list ELSE list FI
			{ $$ = make_if_command ($2, $4, $6); }
	|	IF list THEN list elif_clause FI
			{ $$ = make_if_command ($2, $4, $5); }
	;


group_command:	'{' list '}'
			{ $$ = make_group_command ($2); }
	;

elif_clause:	ELIF list THEN list
			{ $$ = make_if_command ($2, $4, (COMMAND *)NULL); }
	|	ELIF list THEN list ELSE list
			{ $$ = make_if_command ($2, $4, $6); }
	|	ELIF list THEN list elif_clause
			{ $$ = make_if_command ($2, $4, $5); }
	;

case_clause_1:	pattern_list_1
	|	case_clause_sequence pattern_list_1
			{ $2->next = $1; $$ = $2; }
	;

pattern_list_1:	newlines pattern ')' list
			{ $$ = make_pattern_list ($2, $4); }
	|	newlines pattern ')' newlines
			{ $$ = make_pattern_list ($2, (COMMAND *)NULL); }
	|	newlines '(' pattern ')' list
			{ $$ = make_pattern_list ($3, $5); }
	|	newlines '(' pattern ')' newlines
			{ $$ = make_pattern_list ($3, (COMMAND *)NULL); }
	;

case_clause_sequence:  pattern_list
	|	case_clause_sequence pattern_list
			{ $2->next = $1; $$ = $2; }
	;

pattern_list:	newlines pattern ')' list SEMI_SEMI
			{ $$ = make_pattern_list ($2, $4); }
	|	newlines pattern ')' newlines SEMI_SEMI
			{ $$ = make_pattern_list ($2, (COMMAND *)NULL); }
	|	newlines '(' pattern ')' list SEMI_SEMI
			{ $$ = make_pattern_list ($3, $5); }
	|	newlines '(' pattern ')' newlines SEMI_SEMI
			{ $$ = make_pattern_list ($3, (COMMAND *)NULL); }
	;

pattern:	WORD
			{ $$ = make_word_list ($1, (WORD_LIST *)NULL); }
	|	pattern '|' WORD
			{ $$ = make_word_list ($3, $1); }
	;

/* A list allows leading or trailing newlines and
   newlines as operators (equivalent to semicolons).
   It must end with a newline or semicolon.
   Lists are used within commands such as if, for, while.  */

list:		newlines list0
			{
			  $$ = $2;
			  if (need_here_doc)
			    gather_here_documents ();
			 }
	;

list0:		list1
	|	list1 '\n' newlines
	|	list1 '&' newlines
			{
			  if ($1->type == cm_connection)
			    $$ = connect_async_list ($1, (COMMAND *)NULL, '&');
			  else
			    $$ = command_connect ($1, (COMMAND *)NULL, '&');
			}
	|	list1 ';' newlines

	;

list1:		list1 AND_AND newlines list1
			{ $$ = command_connect ($1, $4, AND_AND); }
	|	list1 OR_OR newlines list1
			{ $$ = command_connect ($1, $4, OR_OR); }
	|	list1 '&' newlines list1
			{
			  if ($1->type == cm_connection)
			    $$ = connect_async_list ($1, $4, '&');
			  else
			    $$ = command_connect ($1, $4, '&');
			}
	|	list1 ';' newlines list1
			{ $$ = command_connect ($1, $4, ';'); }
	|	list1 '\n' newlines list1
			{ $$ = command_connect ($1, $4, ';'); }
	|	pipeline
			{ $$ = $1; }
	|	BANG pipeline
			{
			  $2->flags |= CMD_INVERT_RETURN;
			  $$ = $2;
			}
	;

list_terminator:'\n'
	|	';'
	|	yacc_EOF
	;

newlines:
	|	newlines '\n'
	;

/* A simple_list is a list that contains no significant newlines
   and no leading or trailing newlines.  Newlines are allowed
   only following operators, where they are not significant.

   This is what an inputunit consists of.  */

simple_list:	simple_list1
			{
			  $$ = $1;
			  if (need_here_doc)
			    gather_here_documents ();
			}
	|	simple_list1 '&'
			{
			  if ($1->type == cm_connection)
			    $$ = connect_async_list ($1, (COMMAND *)NULL, '&');
			  else
			    $$ = command_connect ($1, (COMMAND *)NULL, '&');
			  if (need_here_doc)
			    gather_here_documents ();
			}
	|	simple_list1 ';'
			{
			  $$ = $1;
			  if (need_here_doc)
			    gather_here_documents ();
			}
	;

simple_list1:	simple_list1 AND_AND newlines simple_list1
			{ $$ = command_connect ($1, $4, AND_AND); }
	|	simple_list1 OR_OR newlines simple_list1
			{ $$ = command_connect ($1, $4, OR_OR); }
	|	simple_list1 '&' simple_list1
			{
			  if ($1->type == cm_connection)
			    $$ = connect_async_list ($1, $3, '&');
			  else
			    $$ = command_connect ($1, $3, '&');
			}
	|	simple_list1 ';' simple_list1
			{ $$ = command_connect ($1, $3, ';'); }
	|	pipeline
			{ $$ = $1; }
	|	BANG pipeline
			{
			  $2->flags |= CMD_INVERT_RETURN;
			  $$ = $2;
			}
	;

pipeline:
		pipeline '|' newlines pipeline
			{ $$ = command_connect ($1, $4, '|'); }
	|	command
			{ $$ = $1; }
	;
%%

/* Initial size to allocate for tokens, and the
   amount to grow them by. */
#define TOKEN_DEFAULT_GROW_SIZE 512

/* The token currently being read. */
static int current_token = 0;

/* The last read token, or NULL.  read_token () uses this for context
   checking. */
static int last_read_token = 0;

/* The token read prior to last_read_token. */
static int token_before_that = 0;

/* If non-zero, it is the token that we want read_token to return
   regardless of what text is (or isn't) present to be read.  This
   is reset by read_token. */
static int token_to_read = 0;

/* Global var is non-zero when end of file has been reached. */
int EOF_Reached = 0;

/* yy_getc () returns the next available character from input or EOF.
   yy_ungetc (c) makes `c' the next character to read.
   init_yy_io (get, unget, type, location) makes the function GET the
   installed function for getting the next character, makes UNGET the
   installed function for un-getting a character, sets the type of stream
   (either string or file) from TYPE, and makes LOCATION point to where
   the input is coming from. */

/* Unconditionally returns end-of-file. */
return_EOF ()
{
  return (EOF);
}

/* Variable containing the current get and unget functions.
   See ./input.h for a clearer description. */
BASH_INPUT bash_input;

/* Set all of the fields in BASH_INPUT to NULL. */
void
initialize_bash_input ()
{
  bash_input.type = 0;
  bash_input.name = (char *)NULL;
  bash_input.location.file = (FILE *)NULL;
  bash_input.location.string = (char *)NULL;
  bash_input.getter = (Function *)NULL;
  bash_input.ungetter = (Function *)NULL;
}

/* Set the contents of the current bash input stream from
   GET, UNGET, TYPE, NAME, and LOCATION. */
void
init_yy_io (get, unget, type, name, location)
     Function *get, *unget;
     int type;
     char *name;
     INPUT_STREAM location;
{
  bash_input.type = type;
  FREE (bash_input.name);

  if (name)
    bash_input.name = savestring (name);
  else
    bash_input.name = (char *)NULL;

#if defined (CRAY)
  memcpy((char *)&bash_input.location.string, (char *)&location.string, sizeof(location));
#else
  bash_input.location = location;
#endif
  bash_input.getter = get;
  bash_input.ungetter = unget;
}

/* Call this to get the next character of input. */
yy_getc ()
{
  return (*(bash_input.getter)) ();
}

/* Call this to unget C.  That is, to make C the next character
   to be read. */
yy_ungetc (c)
     int c;
{
  return (*(bash_input.ungetter)) (c);
}

#if defined (BUFFERED_INPUT)
int
input_file_descriptor ()
{
  switch (bash_input.type)
    {
    case st_stream:
      return (fileno (bash_input.location.file));
    case st_bstream:
      return (bash_input.location.buffered_fd);
    default:
      return (fileno (stdin));
    }
}
#endif /* BUFFERED_INPUT */

/* **************************************************************** */
/*								    */
/*		  Let input be read from readline ().		    */
/*								    */
/* **************************************************************** */

#if defined (READLINE)
char *current_readline_prompt = (char *)NULL;
char *current_readline_line = (char *)NULL;
int current_readline_line_index = 0;

static int
yy_readline_get ()
{
  if (!current_readline_line)
    {
      SigHandler *old_sigint;
      int line_len;

      if (!bash_readline_initialized)
	initialize_readline ();

#if defined (JOB_CONTROL)
      if (job_control)
	give_terminal_to (shell_pgrp);
#endif /* JOB_CONTROL */

      if (signal_is_ignored (SIGINT) == 0)
	{
	  old_sigint = (SigHandler *)set_signal_handler (SIGINT, sigint_sighandler);
	  interrupt_immediately++;
	}

      if (!current_readline_prompt)
	current_readline_line = readline ("");
      else
	current_readline_line = readline (current_readline_prompt);

      if (signal_is_ignored (SIGINT) == 0)
	{
	  interrupt_immediately--;
	  set_signal_handler (SIGINT, old_sigint);
	}

      /* Reset the prompt to whatever is in the decoded value of
	 prompt_string_pointer. */
      reset_readline_prompt ();

      current_readline_line_index = 0;

      if (!current_readline_line)
	return (EOF);

      line_len = strlen (current_readline_line);
      current_readline_line = xrealloc (current_readline_line, 2 + line_len);
      current_readline_line[line_len++] = '\n';
      current_readline_line[line_len] = '\0';
    }

  if (!current_readline_line[current_readline_line_index])
    {
      free (current_readline_line);
      current_readline_line = (char *)NULL;
      return (yy_readline_get ());
    }
  else
    {
      int c = (unsigned char)current_readline_line[current_readline_line_index++];
      return (c);
    }
}

static int
yy_readline_unget (c)
{
  if (current_readline_line_index && current_readline_line)
    current_readline_line[--current_readline_line_index] = c;
  return (c);
}

void  
with_input_from_stdin ()
{
  INPUT_STREAM location;

  if (bash_input.type != st_stdin && stream_on_stack (st_stdin) == 0)
    {
      location.string = current_readline_line;
      init_yy_io (yy_readline_get, yy_readline_unget,
		  st_stdin, "readline stdin", location);
    }
}

#else  /* !READLINE */

void
with_input_from_stdin ()
{
  with_input_from_stream (stdin, "stdin");
}
#endif	/* !READLINE */

/* **************************************************************** */
/*								    */
/*   Let input come from STRING.  STRING is zero terminated.	    */
/*								    */
/* **************************************************************** */

static int
yy_string_get ()
{
  register unsigned char *string;
  register int c;

  string = bash_input.location.string;
  c = EOF;

  /* If the string doesn't exist, or is empty, EOF found. */
  if (string && *string)
    {
      c = *string++;
      bash_input.location.string = string;
    }
  return (c);
}

static int
yy_string_unget (c)
     int c;
{
  *(--bash_input.location.string) = c;
  return (c);
}

void
with_input_from_string (string, name)
     char *string;
     char *name;
{
  INPUT_STREAM location;

  location.string = string;

  init_yy_io (yy_string_get, yy_string_unget, st_string, name, location);
}

/* **************************************************************** */
/*								    */
/*		     Let input come from STREAM.		    */
/*								    */
/* **************************************************************** */

static int
yy_stream_get ()
{
  int result = EOF;

  if (bash_input.location.file)
#if defined (NO_READ_RESTART_ON_SIGNAL)
    result = (unsigned char)getc_with_restart (bash_input.location.file);
#else
    result = (unsigned char)getc (bash_input.location.file);
#endif /* !NO_READ_RESTART_ON_SIGNAL */
  return (result);
}

static int
yy_stream_unget (c)
     int c;
{
#if defined (NO_READ_RESTART_ON_SIGNAL)
  return (ungetc_with_restart (c, bash_input.location.file));
#else
  return (ungetc (c, bash_input.location.file));
#endif
}

void
with_input_from_stream (stream, name)
     FILE *stream;
     char *name;
{
  INPUT_STREAM location;

  location.file = stream;
  init_yy_io (yy_stream_get, yy_stream_unget, st_stream, name, location);
}

typedef struct stream_saver {
  struct stream_saver *next;
  BASH_INPUT bash_input;
  int line;
#if defined (BUFFERED_INPUT)
  BUFFERED_STREAM *bstream;
#endif /* BUFFERED_INPUT */
} STREAM_SAVER;

/* The globally known line number. */
int line_number = 0;

STREAM_SAVER *stream_list = (STREAM_SAVER *)NULL;

push_stream ()
{
  STREAM_SAVER *saver = (STREAM_SAVER *)xmalloc (sizeof (STREAM_SAVER));

  xbcopy ((char *)&bash_input, (char *)&(saver->bash_input), sizeof (BASH_INPUT));

#if defined (BUFFERED_INPUT)
  saver->bstream = (BUFFERED_STREAM *)NULL;
  /* If we have a buffered stream, clear out buffers[fd]. */
  if (bash_input.type == st_bstream && bash_input.location.buffered_fd >= 0)
    {
      saver->bstream = buffers[bash_input.location.buffered_fd];
      buffers[bash_input.location.buffered_fd] = (BUFFERED_STREAM *)NULL;
    }
#endif /* BUFFERED_INPUT */

  saver->line = line_number;
  bash_input.name = (char *)NULL;
  saver->next = stream_list;
  stream_list = saver;
  EOF_Reached = line_number = 0;
}

pop_stream ()
{
  int temp;

  if (!stream_list)
    EOF_Reached = 1;
  else
    {
      STREAM_SAVER *saver = stream_list;

      EOF_Reached = 0;
      stream_list = stream_list->next;

      init_yy_io (saver->bash_input.getter,
		  saver->bash_input.ungetter,
		  saver->bash_input.type,
		  saver->bash_input.name,
		  saver->bash_input.location);

#if defined (BUFFERED_INPUT)
      /* If we have a buffered stream, restore buffers[fd]. */
      /* If the input file descriptor was changed while this was on the
	 save stack, update the buffered fd to the new file descriptor and
	 re-establish the buffer <-> bash_input fd correspondence. */
      if (bash_input.type == st_bstream && bash_input.location.buffered_fd >= 0)
        {
          if (bash_input_fd_changed)
	    {
	      bash_input_fd_changed = 0;
	      if (default_buffered_input >= 0)
		{
		  bash_input.location.buffered_fd = default_buffered_input;
		  saver->bstream->b_fd = default_buffered_input;
		}
	    }
	  buffers[bash_input.location.buffered_fd] = saver->bstream;
        }
#endif /* BUFFERED_INPUT */

      line_number = saver->line;

      FREE (saver->bash_input.name);
      free (saver);
    }
}

/* Return 1 if a stream of type TYPE is saved on the stack. */
int
stream_on_stack (type)
     int type;
{
  register STREAM_SAVER *s;
 
  for (s = stream_list; s; s = s->next)
    if (s->bash_input.type == type)
      return 1;
  return 0;
}


/*
 * This is used to inhibit alias expansion and reserved word recognition
 * inside case statement pattern lists.  A `case statement pattern list'
 * is:
 *	everything between the `in' in a `case word in' and the next ')'
 *	or `esac'
 *	everything between a `;;' and the next `)' or `esac'
 */
static int in_case_pattern_list = 0;

#if defined (ALIAS)
/*
 * Pseudo-global variables used in implementing token-wise alias expansion.
 */

static int expand_next_token = 0;

/*
 * Pushing and popping strings.  This works together with shell_getc to 
 * implement alias expansion on a per-token basis.
 */

typedef struct string_saver {
  struct string_saver *next;
  int expand_alias;  /* Value to set expand_alias to when string is popped. */
  char *saved_line;
  int saved_line_size, saved_line_index, saved_line_terminator;
} STRING_SAVER;

STRING_SAVER *pushed_string_list = (STRING_SAVER *)NULL;

static void save_expansion ();

/*
 * Push the current shell_input_line onto a stack of such lines and make S
 * the current input.  Used when expanding aliases.  EXPAND is used to set
 * the value of expand_next_token when the string is popped, so that the
 * word after the alias in the original line is handled correctly when the
 * alias expands to multiple words.  TOKEN is the token that was expanded
 * into S; it is saved and used to prevent infinite recursive expansion.
 */
static void
push_string (s, expand, token)
     char *s;
     int expand;
     char *token;
{
  STRING_SAVER *temp = (STRING_SAVER *) xmalloc (sizeof (STRING_SAVER));

  temp->expand_alias = expand;
  temp->saved_line = shell_input_line;
  temp->saved_line_size = shell_input_line_size;
  temp->saved_line_index = shell_input_line_index;
  temp->saved_line_terminator = shell_input_line_terminator;
  temp->next = pushed_string_list;
  pushed_string_list = temp;

  save_expansion (token);

  shell_input_line = s;
  shell_input_line_size = strlen (s);
  shell_input_line_index = 0;
  shell_input_line_terminator = '\0';
  expand_next_token = 0;
}

/*
 * Make the top of the pushed_string stack be the current shell input.
 * Only called when there is something on the stack.  Called from shell_getc
 * when it thinks it has consumed the string generated by an alias expansion
 * and needs to return to the original input line.
 */
static void
pop_string ()
{
  STRING_SAVER *t;

  FREE (shell_input_line);
  shell_input_line = pushed_string_list->saved_line;
  shell_input_line_index = pushed_string_list->saved_line_index;
  shell_input_line_size = pushed_string_list->saved_line_size;
  shell_input_line_terminator = pushed_string_list->saved_line_terminator;
  expand_next_token = pushed_string_list->expand_alias;

  t = pushed_string_list;
  pushed_string_list = pushed_string_list->next;
  free((char *)t);
}

static void
free_string_list ()
{
  register STRING_SAVER *t = pushed_string_list, *t1;

  while (t)
    {
      t1 = t->next;
      FREE (t->saved_line);
      free ((char *)t);
      t = t1;
    }
  pushed_string_list = (STRING_SAVER *)NULL;
}

/* This is a stack to save the values of all tokens for which alias
   expansion has been performed during the current call to read_token ().
   It is used to prevent alias expansion loops:

      alias foo=bar
      alias bar=baz
      alias baz=foo

   Ideally this would be taken care of by push and pop string, but because
   of when strings are popped the stack will not contain the correct
   strings to test against.  (The popping is done in shell_getc, so that when
   the current string is exhausted, shell_getc can simply pop that string off
   the stack, restore the previous string, and continue with the character
   following the token whose expansion was originally pushed on the stack.)

   What we really want is a record of all tokens that have been expanded for
   aliases during the `current' call to read_token().  This does that, at the
   cost of being somewhat special-purpose (OK, OK vile and unclean). */

typedef struct _exp_saver {
      struct _exp_saver *next;
      char *saved_token;
} EXPANSION_SAVER;

EXPANSION_SAVER *expanded_token_stack = (EXPANSION_SAVER *)NULL;

static void
save_expansion (s)
     char *s;
{
  EXPANSION_SAVER *t;

  t = (EXPANSION_SAVER *) xmalloc (sizeof (EXPANSION_SAVER));
  t->saved_token = savestring (s);
  t->next = expanded_token_stack;
  expanded_token_stack = t;
}

/* Return 1 if TOKEN has already been expanded in the current `stack' of
   expansions.  If it has been expanded already, it will appear as the value
   of saved_token for some entry in the stack of expansions created for the
   current token being expanded. */
static int
token_has_been_expanded (token)
     char *token;
{
  register EXPANSION_SAVER *t = expanded_token_stack;

  while (t)
    {
      if (STREQ (token, t->saved_token))
	return (1);
      t = t->next;
    }
  return (0);
}

static void
free_expansion_stack ()
{
  register EXPANSION_SAVER *t = expanded_token_stack, *t1;

  while (t)
    {
      t1 = t->next;
      free (t->saved_token);
      free (t);
      t = t1;
    }
  expanded_token_stack = (EXPANSION_SAVER *)NULL;
}

#endif /* ALIAS */

/* Return a line of text, taken from wherever yylex () reads input.
   If there is no more input, then we return NULL.  If REMOVE_QUOTED_NEWLINE
   is non-zero, we remove unquoted \<newline> pairs.  This is used by
   read_secondary_line to read here documents. */
static char *
read_a_line (remove_quoted_newline)
     int remove_quoted_newline;
{
  static char *line_buffer = (char *)NULL;
  static int buffer_size = 0;
  int indx = 0, c, peekc, pass_next;

  pass_next = 0;
  while (1)
    {
      c = yy_getc ();

      /* Allow immediate exit if interrupted during input. */
      QUIT;

      if (c == 0)
	continue;

      /* If there is no more input, then we return NULL. */
      if (c == EOF)
	{
	  if (indx == 0)
	    return ((char *)NULL);
	  c = '\n';
	}

      /* `+2' in case the final character in the buffer is a newline. */
      if (indx + 2 > buffer_size)
	if (!buffer_size)
	  line_buffer = xmalloc (buffer_size = 128);
	else
	  line_buffer = xrealloc (line_buffer, buffer_size += 128);

      /* IF REMOVE_QUOTED_NEWLINES is non-zero, we are reading a
	 here document with an unquoted delimiter.  In this case,
	 the line will be expanded as if it were in double quotes.
	 We allow a backslash to escape the next character, but we
	 need to treat the backslash specially only if a backslash
	 quoting a backslash-newline pair appears in the line. */
      if (pass_next)
        {
	  line_buffer[indx++] = c;
	  pass_next = 0;
        }
      else if (c == '\\' && remove_quoted_newline)
	{
	  peekc = yy_getc ();
	  if (peekc == '\n')
	    continue;	/* Make the unquoted \<newline> pair disappear. */
	  else
	    {
	      yy_ungetc (peekc);
	      pass_next = 1;
	      line_buffer[indx++] = c;		/* Preserve the backslash. */
	    }
	}
      else
	line_buffer[indx++] = c;

      if (c == '\n')
	{
	  line_buffer[indx] = '\0';
	  return (line_buffer);
	}
    }
}

/* Return a line as in read_a_line (), but insure that the prompt is
   the secondary prompt.  This is used to read the lines of a here
   document.  REMOVE_QUOTED_NEWLINE is non-zero if we should remove
   newlines quoted with backslashes while reading the line.  It is
   non-zero unless the delimiter of the here document was quoted. */
char *
read_secondary_line (remove_quoted_newline)
     int remove_quoted_newline;
{
  prompt_string_pointer = &ps2_prompt;
  prompt_again ();
  return (read_a_line (remove_quoted_newline));
}


/* **************************************************************** */
/*								    */
/*				YYLEX ()			    */
/*								    */
/* **************************************************************** */

/* Reserved words.  These are only recognized as the first word of a
   command. */
STRING_INT_ALIST word_token_alist[] = {
  { "if", IF },
  { "then", THEN },
  { "else", ELSE },
  { "elif", ELIF },
  { "fi", FI },
  { "case", CASE },
  { "esac", ESAC },
  { "for", FOR },
#if defined (SELECT_COMMAND)
  { "select", SELECT },
#endif
  { "while", WHILE },
  { "until", UNTIL },
  { "do", DO },
  { "done", DONE },
  { "in", IN },
  { "function", FUNCTION },
  { "{", '{' },
  { "}", '}' },
  { "!", BANG },
  { (char *)NULL, 0}
};

/* Return the next shell input character.  This always reads characters
   from shell_input_line; when that line is exhausted, it is time to
   read the next line.  This is called by read_token when the shell is
   processing normal command input. */
static int
shell_getc (remove_quoted_newline)
     int remove_quoted_newline;
{
  int c;

  QUIT;

#if defined (ALIAS)
  /* If shell_input_line[shell_input_line_index] == 0, but there is
     something on the pushed list of strings, then we don't want to go
     off and get another line.  We let the code down below handle it. */

  if (!shell_input_line || ((!shell_input_line[shell_input_line_index]) &&
			    (pushed_string_list == (STRING_SAVER *)NULL)))
#else /* !ALIAS */
  if (!shell_input_line || !shell_input_line[shell_input_line_index])
#endif /* !ALIAS */
    {
      register int i, l;

      restart_read_next_line:

      line_number++;

    restart_read:

      /* Allow immediate exit if interrupted during input. */
      QUIT;

      i = 0;
      shell_input_line_terminator = 0;

#if defined (JOB_CONTROL)
      /* This can cause a problem when reading a command as the result
	 of a trap, when the trap is called from flush_child.  This call
	 had better not cause jobs to disappear from the job table in
	 that case, or we will have big trouble. */
      notify_and_cleanup ();
#else /* !JOB_CONTROL */
      cleanup_dead_jobs ();
#endif /* !JOB_CONTROL */

#if defined (READLINE)
      if (interactive && bash_input.type != st_string && no_line_editing)
#else
      if (interactive && bash_input.type != st_string)
#endif
	print_prompt ();

      if (bash_input.type == st_stream)
	clearerr (stdin);

      while (c = yy_getc ())
	{
	  /* Allow immediate exit if interrupted during input. */
	  QUIT;

	  if (i + 2 > shell_input_line_size)
	    shell_input_line =
	      xrealloc (shell_input_line, shell_input_line_size += 256);

	  if (c == EOF)
	    {
	      if (bash_input.type == st_stream)
		clearerr (stdin);

	      if (!i)
		shell_input_line_terminator = EOF;

	      shell_input_line[i] = '\0';
	      break;
	    }

	  shell_input_line[i++] = c;

	  if (c == '\n')
	    {
	      shell_input_line[--i] = '\0';
	      current_command_line_count++;
	      break;
	    }
	}
      shell_input_line_index = 0;
      shell_input_line_len = i;		/* == strlen (shell_input_line) */

#if defined (HISTORY)
      if (interactive && shell_input_line && shell_input_line[0])
	{
	  char *expansions;

	  expansions = pre_process_line (shell_input_line, 1, 1);

	  free (shell_input_line);
	  shell_input_line = expansions;
	  shell_input_line_len = shell_input_line ?
				 strlen (shell_input_line) :
				 0;
	  if (!shell_input_line_len)
	    current_command_line_count--;

	  /* We have to force the xrealloc below because we don't know the
	     true allocated size of shell_input_line anymore. */
	  shell_input_line_size = shell_input_line_len;
	}
#endif /* HISTORY */

      if (shell_input_line)
	{
	  /* Lines that signify the end of the shell's input should not be
	     echoed. */
	  if (echo_input_at_read && (shell_input_line[0] ||
				     shell_input_line_terminator != EOF))
	    fprintf (stderr, "%s\n", shell_input_line);
	}
      else
	{
	  shell_input_line_size = 0;
	  prompt_string_pointer = &current_prompt_string;
	  prompt_again ();
	  goto restart_read;
	}

      /* Add the newline to the end of this string, iff the string does
	 not already end in an EOF character.  */
      if (shell_input_line_terminator != EOF)
	{
	  l = shell_input_line_len;	/* was a call to strlen */

	  if (l + 3 > shell_input_line_size)
	    shell_input_line = xrealloc (shell_input_line,
					1 + (shell_input_line_size += 2));

	  shell_input_line[l] = '\n';
	  shell_input_line[l + 1] = '\0';
	}
    }
  
  c = shell_input_line[shell_input_line_index];

  if (c)
    shell_input_line_index++;

  if (c == '\\' && remove_quoted_newline &&
      shell_input_line[shell_input_line_index] == '\n')
    {
	prompt_again ();
	goto restart_read_next_line;
    }

#if defined (ALIAS)
  /* If C is NULL, we have reached the end of the current input string.  If
     pushed_string_list is non-empty, it's time to pop to the previous string
     because we have fully consumed the result of the last alias expansion.
     Do it transparently; just return the next character of the string popped
     to. */
  if (!c && (pushed_string_list != (STRING_SAVER *)NULL))
    {
      pop_string ();
      c = shell_input_line[shell_input_line_index];
      if (c)
	shell_input_line_index++;
    }
#endif /* ALIAS */

  if (!c && shell_input_line_terminator == EOF)
    {
      if (shell_input_line_index != 0)
	return ('\n');
      else
	return (EOF);
    }

  return ((unsigned char)c);
}

/* Put C back into the input for the shell. */
static void
shell_ungetc (c)
     int c;
{
  if (shell_input_line && shell_input_line_index)
    shell_input_line[--shell_input_line_index] = c;
}

/* Discard input until CHARACTER is seen. */
static void
discard_until (character)
     int character;
{
  int c;

  while ((c = shell_getc (0)) != EOF && c != character)
    ;

  if (c != EOF)
    shell_ungetc (c);
}

/* Place to remember the token.  We try to keep the buffer
   at a reasonable size, but it can grow. */
static char *token = (char *)NULL;

/* Current size of the token buffer. */
static int token_buffer_size = 0;

void
execute_prompt_command (command)
     char *command;
{
  Function *temp_last, *temp_this;
  char *last_lastarg;
  int temp_exit_value, temp_eof_encountered;

  temp_last = last_shell_builtin;
  temp_this = this_shell_builtin;
  temp_exit_value = last_command_exit_value;
  temp_eof_encountered = eof_encountered;
  last_lastarg = get_string_value ("_");
  if (last_lastarg)
    last_lastarg = savestring (last_lastarg);

  parse_and_execute (savestring (command), "PROMPT_COMMAND", 0);

  last_shell_builtin = temp_last;
  this_shell_builtin = temp_this;
  last_command_exit_value = temp_exit_value;
  eof_encountered = temp_eof_encountered;

  bind_variable ("_", last_lastarg);
  FREE (last_lastarg);

  if (token_to_read == '\n')
    token_to_read = 0;
}

/* Command to read_token () explaining what we want it to do. */
#define READ 0
#define RESET 1
#define prompt_is_ps1 \
      (!prompt_string_pointer || prompt_string_pointer == &ps1_prompt)

/* Function for yyparse to call.  yylex keeps track of
   the last two tokens read, and calls read_token.  */

yylex ()
{
  if (interactive && (!current_token || current_token == '\n'))
    {
      /* Before we print a prompt, we might have to check mailboxes.
	 We do this only if it is time to do so. Notice that only here
	 is the mail alarm reset; nothing takes place in check_mail ()
	 except the checking of mail.  Please don't change this. */
      if (prompt_is_ps1 && time_to_check_mail ())
	{
	  check_mail ();
	  reset_mail_timer ();
	}

      /* Avoid printing a prompt if we're not going to read anything, e.g.
	 after resetting the parser with read_token (RESET). */
      if (token_to_read == 0 && interactive)
	prompt_again ();
    }

  token_before_that = last_read_token;
  last_read_token = current_token;
  current_token = read_token (READ);
  return (current_token);
}

/* Called from shell.c when Control-C is typed at top level.  Or
   by the error rule at top level. */
reset_parser ()
{
  read_token (RESET);
}
  
/* When non-zero, we have read the required tokens
   which allow ESAC to be the next one read. */
static int allow_esac_as_next = 0;

/* When non-zero, accept single '{' as a token itself. */
static int allow_open_brace = 0;

/* DELIMITERS is a stack of the nested delimiters that we have
   encountered so far. */
static char *delimiters = (char *)NULL;

/* Offset into the stack of delimiters. */
int delimiter_depth = 0;

/* How many slots are allocated to DELIMITERS. */
static int delimiter_space = 0;

void
gather_here_documents ()
{
  int r = 0;
  while (need_here_doc)
    {
      make_here_document (redir_stack[r++]);
      need_here_doc--;
    }
}

/* Macro for accessing the top delimiter on the stack.  Returns the
   delimiter or zero if none. */
#define current_delimiter() \
  (delimiter_depth ? delimiters[delimiter_depth - 1] : 0)

#define push_delimiter(character) \
  do \
    { \
      if (delimiter_depth + 2 > delimiter_space) \
	delimiters = xrealloc \
	  (delimiters, (delimiter_space += 10) * sizeof (char)); \
      delimiters[delimiter_depth] = character; \
      delimiter_depth++; \
    } \
  while (0)

/* When non-zero, an open-brace used to create a group is awaiting a close
   brace partner. */
static int open_brace_awaiting_satisfaction = 0;

#define command_token_position(token) \
  (((token) == ASSIGNMENT_WORD) || \
   ((token) != SEMI_SEMI && reserved_word_acceptable(token)))

#define assignment_acceptable(token) command_token_position(token) && \
					(in_case_pattern_list == 0)

/* Check to see if TOKEN is a reserved word and return the token
   value if it is. */
#define CHECK_FOR_RESERVED_WORD(tok) \
  do { \
    if (!dollar_present && !quoted && \
	reserved_word_acceptable (last_read_token)) \
      { \
	int i; \
	for (i = 0; word_token_alist[i].word != (char *)NULL; i++) \
	  if (STREQ (tok, word_token_alist[i].word)) \
	    { \
	      if (in_case_pattern_list && (word_token_alist[i].token != ESAC)) \
		break; \
\
	      if (word_token_alist[i].token == ESAC) \
		in_case_pattern_list = 0; \
\
	      if (word_token_alist[i].token == '{') \
		open_brace_awaiting_satisfaction++; \
\
	      if (word_token_alist[i].token == '}' && open_brace_awaiting_satisfaction) \
		open_brace_awaiting_satisfaction--; \
\
	      return (word_token_alist[i].token); \
	    } \
      } \
  } while (0)

/* Read the next token.  Command can be READ (normal operation) or 
   RESET (to normalize state). */
static int
read_token (command)
     int command;
{
  int character;		/* Current character. */
  int peek_char;		/* Temporary look-ahead character. */
  int result;			/* The thing to return. */
  WORD_DESC *the_word;		/* The value for YYLVAL when a WORD is read. */

  if (token_buffer_size < TOKEN_DEFAULT_GROW_SIZE)
    {
      FREE (token);
      token = xmalloc (token_buffer_size = TOKEN_DEFAULT_GROW_SIZE);
    }

  if (command == RESET)
    {
      delimiter_depth = 0;	/* No delimiters found so far. */
      open_brace_awaiting_satisfaction = 0;
      in_case_pattern_list = 0;

#if defined (ALIAS)
      if (pushed_string_list)
	{
	  free_string_list ();
	  pushed_string_list = (STRING_SAVER *)NULL;
	}

      if (expanded_token_stack)
	{
	  free_expansion_stack ();
	  expanded_token_stack = (EXPANSION_SAVER *)NULL;
	}

      expand_next_token = 0;
#endif /* ALIAS */

      if (shell_input_line)
	{
	  free (shell_input_line);
	  shell_input_line = (char *)NULL;
	  shell_input_line_size = shell_input_line_index = 0;
	}
      last_read_token = '\n';
      token_to_read = '\n';
      return ('\n');
    }

  if (token_to_read)
    {
      int rt = token_to_read;
      token_to_read = 0;
      return (rt);
    }

#if defined (ALIAS)
  /* If we hit read_token () and there are no saved strings on the
     pushed_string_list, then we are no longer currently expanding a
     token.  This can't be done in pop_stream, because pop_stream
     may pop the stream before the current token has finished being
     completely expanded (consider what happens when we alias foo to foo,
     and then try to expand it). */
  if (!pushed_string_list && expanded_token_stack)
    {
      free_expansion_stack ();
      expanded_token_stack = (EXPANSION_SAVER *)NULL;
    }

  /* This is a place to jump back to once we have successfully expanded a
     token with an alias and pushed the string with push_string () */
 re_read_token:

#endif /* ALIAS */

  /* Read a single word from input.  Start by skipping blanks. */
  while ((character = shell_getc (1)) != EOF && whitespace (character));

  if (character == EOF)
    {
      EOF_Reached = 1;
      return (yacc_EOF);
    }

  if (character == '#' && (!interactive || interactive_comments))
    {
      /* A comment.  Discard until EOL or EOF, and then return a newline. */
      discard_until ('\n');
      shell_getc (0);

      /* If we're about to return an unquoted newline, we can go and collect
	 the text of any pending here documents. */
      if (need_here_doc)
        gather_here_documents ();

#if defined (ALIAS)
      expand_next_token = 0;
#endif /* ALIAS */

      return ('\n');
    }

  if (character == '\n')
    {
      /* If we're about to return an unquoted newline, we can go and collect
	 the text of any pending here document. */
      if (need_here_doc)
	gather_here_documents ();

#if defined (ALIAS)
      expand_next_token = 0;
#endif /* ALIAS */

      return (character);
    }

  if (member (character, "()<>;&|"))
    {
#if defined (ALIAS)
      /* Turn off alias tokenization iff this character sequence would
	 not leave us ready to read a command. */
      if (character == '<' || character == '>')
	expand_next_token = 0;
#endif /* ALIAS */

      /* Please note that the shell does not allow whitespace to
	 appear in between tokens which are character pairs, such as
	 "<<" or ">>".  I believe this is the correct behaviour. */
      if (character == (peek_char = shell_getc (1)))
	{
	  switch (character)
	    {
	      /* If '<' then we could be at "<<" or at "<<-".  We have to
		 look ahead one more character. */
	    case '<':
	      peek_char = shell_getc (1);
	      if (peek_char == '-')
		return (LESS_LESS_MINUS);
	      else
		{
		  shell_ungetc (peek_char);
		  return (LESS_LESS);
		}

	    case '>':
	      return (GREATER_GREATER);

	    case ';':
	      in_case_pattern_list = 1;
#if defined (ALIAS)
	      expand_next_token = 0;
#endif /* ALIAS */
	      return (SEMI_SEMI);

	    case '&':
	      return (AND_AND);

	    case '|':
	      return (OR_OR);
	    }
	}
      else
	{
	  if (peek_char == '&')
	    {
	      switch (character)
		{
		case '<': return (LESS_AND);
		case '>': return (GREATER_AND);
		}
	    }
	  if (character == '<' && peek_char == '>')
	    return (LESS_GREATER);
	  if (character == '>' && peek_char == '|')
	    return (GREATER_BAR);
	  if (peek_char == '>' && character == '&')
	    return (AND_GREATER);
	}
      shell_ungetc (peek_char);

      /* If we look like we are reading the start of a function
	 definition, then let the reader know about it so that
	 we will do the right thing with `{'. */
      if (character == ')' &&
	  last_read_token == '(' && token_before_that == WORD)
	{
	  allow_open_brace = 1;
#if defined (ALIAS)
	  expand_next_token = 0;
#endif /* ALIAS */
	}

      if (in_case_pattern_list && (character == ')'))
	in_case_pattern_list = 0;

#if defined (PROCESS_SUBSTITUTION)
      /* Check for the constructs which introduce process substitution.
	 Shells running in `posix mode' don't do process substitution. */
      if (posixly_correct ||
	  (((character == '>' || character == '<') && peek_char == '(') == 0))
#endif /* PROCESS_SUBSTITUTION */
	return (character);
    }

  /* Hack <&- (close stdin) case. */
  if (character == '-')
    {
      switch (last_read_token)
	{
	case LESS_AND:
	case GREATER_AND:
	  return (character);
	}
    }
  
  /* Okay, if we got this far, we have to read a word.  Read one,
     and then check it against the known ones. */
  {
    /* Index into the token that we are building. */
    int token_index = 0;

    /* ALL_DIGITS becomes zero when we see a non-digit. */
    int all_digits = digit (character);

    /* DOLLAR_PRESENT becomes non-zero if we see a `$'. */
    int dollar_present = 0;

    /* QUOTED becomes non-zero if we see one of ("), ('), (`), or (\). */
    int quoted = 0;

    /* Non-zero means to ignore the value of the next character, and just
       to add it no matter what. */
    int pass_next_character = 0;

    /* Non-zero means parsing a dollar-paren construct.  It is the count of
       un-quoted closes we need to see. */
    int dollar_paren_level = 0;

    /* Non-zero means parsing a dollar-bracket construct ($[...]).  It is
       the count of un-quoted `]' characters we need to see. */
    int dollar_bracket_level = 0;

    /* Non-zero means parsing a `${' construct.  It is the count of
       un-quoted `}' we need to see. */
    int dollar_brace_level = 0;

    /* A level variable for parsing '${ ... }' constructs inside of double
       quotes. */
    int delimited_brace_level = 0;

    /* A boolean variable denoting whether or not we are currently parsing
       a double-quoted string embedded in a $( ) or ${ } construct. */
    int embedded_quoted_string = 0;

    /* Another level variable.  This one is for dollar_parens inside of
       double-quotes. */
    int delimited_paren_level = 0;

    /* The current delimiting character. */
    int cd;

    for (;;)
      {
	if (character == EOF)
	  goto got_token;

	if (pass_next_character)
	  {
	    pass_next_character = 0;
	    goto got_character;
	  }

	cd = current_delimiter ();

	if (cd && character == '\\' && cd != '\'')
	  {
	    peek_char = shell_getc (0);
	    if (peek_char != '\\')
	      shell_ungetc (peek_char);
	    else
	      {
		token[token_index++] = character;
		goto got_character;
	      }
	  }

	/* Handle backslashes.  Quote lots of things when not inside of
	   double-quotes, quote some things inside of double-quotes. */
	   
	if (character == '\\' && (!delimiter_depth || cd != '\''))
	  {
	    peek_char = shell_getc (0);

	    /* Backslash-newline is ignored in all cases excepting
	       when quoted with single quotes. */
	    if (peek_char == '\n')
	      {
		character = '\n';
		goto next_character;
	      }
	    else
	      {
		shell_ungetc (peek_char);

		/* If the next character is to be quoted, do it now. */
		if (!cd || cd == '`' ||
		    (cd == '"' && member (peek_char, slashify_in_quotes)))
		  {
		    pass_next_character++;
		    quoted = 1;
		    goto got_character;
		  }
	      }
	  }

	/* This is a hack, in its present form.  If a backquote substitution
	   appears within double quotes, everything within the backquotes
	   should be read as part of a single word.  Jesus.  Now I see why
	   Korn introduced the $() form. */
	if (delimiter_depth && (cd == '"') && (character == '`'))
	  {
	    push_delimiter (character);
	    goto got_character;
	  }

	cd = current_delimiter ();		/* XXX - may not need */
	if (delimiter_depth)
	  {
	    if (character == cd)
	      {
	      	/* If we see a double quote while parsing a double-quoted
		  $( ) or ${ }, and we have not seen ) or }, respectively,
	      	   note that we are in the middle of reading an embedded
		   quoted string. */
		if ((delimited_paren_level || delimited_brace_level) &&
		    (character == '"'))
		  {
		    embedded_quoted_string = !embedded_quoted_string;
		    goto got_character;
		  }
		
		delimiter_depth--;
		goto got_character;
	      }
	  }

	if (cd != '\'')
	  {
#if defined (PROCESS_SUBSTITUTION)
	    if (character == '$' || character == '<' || character == '>')
#else
	    if (character == '$')
#endif /* !PROCESS_SUBSTITUTION */
	      {
	      	/* If we're in the middle of parsing a $( ) or ${ }
	      	   construct with an embedded quoted string, don't
	      	   bother looking at this character any further. */
	      	if (embedded_quoted_string)
	      	  goto got_character;

		peek_char = shell_getc (1);
		shell_ungetc (peek_char);
		if (peek_char == '(')
		  {
		    if (!delimiter_depth)
		      dollar_paren_level++;
		    else
		      delimited_paren_level++;

		    pass_next_character++;
		    goto got_character;
		  }
		else if (peek_char == '[' && character == '$')
		  {
		    if (!delimiter_depth)
		      dollar_bracket_level++;

		    pass_next_character++;
		    goto got_character;
		  }
		/* This handles ${...} constructs. */
		else if (peek_char == '{' && character == '$')
		  {
		    if (!delimiter_depth)
		      dollar_brace_level++;
		    else
		      delimited_brace_level++;

		    pass_next_character++;
		    goto got_character;
		  }
	      }

	    /* If we are parsing a $() or $[] construct, we need to balance
	       parens and brackets inside the construct.  This whole function
	       could use a rewrite. */
	    if (character == '(' && !embedded_quoted_string)
	      {
		if (delimiter_depth && delimited_paren_level)
		  delimited_paren_level++;

		if (!delimiter_depth && dollar_paren_level)
		  dollar_paren_level++;
	      }

	    if (character == '[')
	      {
		if (!delimiter_depth && dollar_bracket_level)
		  dollar_bracket_level++;
	      }

	    if (character == '{' && !embedded_quoted_string)
	      {
	      	if (delimiter_depth && delimited_brace_level)
	      	  delimited_brace_level++;

	      	if (!delimiter_depth && dollar_brace_level)
	      	  dollar_brace_level++;
	      }

	    /* This code needs to take into account whether we are inside a
	       case statement pattern list, and whether this paren is supposed
	       to terminate it (hey, it could happen).  It's not as simple
	       as just using in_case_pattern_list, because we're not parsing
	       anything while we're reading a $( ) construct.  Maybe we
	       should move that whole mess into the yacc parser. */
	    if (character == ')' && !embedded_quoted_string)
	      {
		if (delimiter_depth && delimited_paren_level)
		  delimited_paren_level--;

		if (!delimiter_depth && dollar_paren_level)
		  {
		    dollar_paren_level--;
		    goto got_character;
		  }
	      }

	    if (character == ']')
	      {
		if (!delimiter_depth && dollar_bracket_level)
		  {
		    dollar_bracket_level--;
		    goto got_character;
		  }
	      }

	    if (character == '}' && !embedded_quoted_string)
	      {
		if (delimiter_depth && delimited_brace_level)
		  delimited_brace_level--;

		if (!delimiter_depth && dollar_brace_level)
		  {
		    dollar_brace_level--;
		    goto got_character;
		  }
	      }
	  }

	if (!dollar_paren_level && !dollar_bracket_level &&
	    !dollar_brace_level && !delimiter_depth &&
	    member (character, " \t\n;&()|<>"))
	  {
	    shell_ungetc (character);
	    goto got_token;
	  }
    
	if (!delimiter_depth)
	  {
	    if (character == '"' || character == '`' || character == '\'')
	      {
		push_delimiter (character);

		quoted = 1;
		goto got_character;
	      }
	  }

	if (all_digits)
	  all_digits = digit (character);
	if (character == '$')
	  dollar_present = 1;

      got_character:

	if (character == CTLESC || character == CTLNUL)
	  token[token_index++] = CTLESC;

	token[token_index++] = character;

	if (token_index == (token_buffer_size - 1))
	  {
	    token_buffer_size += TOKEN_DEFAULT_GROW_SIZE;
	    token = xrealloc (token, token_buffer_size);
	  }
	next_character:
	if (character == '\n' && interactive && bash_input.type != st_string)
	  prompt_again ();

	/* We want to remove quoted newlines (that is, a \<newline> pair)
	   unless we are within single quotes or pass_next_character is
	   set (the shell equivalent of literal-next). */
	character = shell_getc
	  ((current_delimiter () != '\'') && (!pass_next_character));
      }

  got_token:

    token[token_index] = '\0';
	
    if ((delimiter_depth || dollar_paren_level || dollar_bracket_level) &&
	character == EOF)
      {
	char reporter = '\0';

	if (!delimiter_depth)
	  {
	    if (dollar_paren_level)
	      reporter = ')';
	    else if (dollar_bracket_level)
	      reporter = ']';
	  }

	if (!reporter)
	  reporter = current_delimiter ();

	report_error ("unexpected EOF while looking for `%c'", reporter);
	return (-1);
      }

    if (all_digits)
      {
	/* Check to see what thing we should return.  If the last_read_token
	   is a `<', or a `&', or the character which ended this token is
	   a '>' or '<', then, and ONLY then, is this input token a NUMBER.
	   Otherwise, it is just a word, and should be returned as such. */

	if (character == '<' || character == '>' ||
	    last_read_token == LESS_AND || last_read_token == GREATER_AND)
	  {
	    yylval.number = atoi (token);
	    return (NUMBER);
	  }
      }

    /* Handle special case.  IN is recognized if the last token
       was WORD and the token before that was FOR or CASE. */
    if ((last_read_token == WORD) &&
#if defined (SELECT_COMMAND)
	((token_before_that == FOR) || (token_before_that == CASE) || (token_before_that == SELECT)) &&
#else
	((token_before_that == FOR) || (token_before_that == CASE)) &&
#endif
	(token[0] == 'i' && token[1] == 'n' && !token[2]))
      {
	if (token_before_that == CASE)
	  {
	    in_case_pattern_list = 1;
	    allow_esac_as_next++;
	  }
	return (IN);
      }

    /* Ditto for DO in the FOR case. */
#if defined (SELECT_COMMAND)
    if ((last_read_token == WORD) && ((token_before_that == FOR) || (token_before_that == SELECT)) &&
#else
    if ((last_read_token == WORD) && (token_before_that == FOR) &&
#endif
	(token[0] == 'd' && token[1] == 'o' && !token[2]))
      return (DO);

    /* Ditto for ESAC in the CASE case. 
       Specifically, this handles "case word in esac", which is a legal
       construct, certainly because someone will pass an empty arg to the
       case construct, and we don't want it to barf.  Of course, we should
       insist that the case construct has at least one pattern in it, but
       the designers disagree. */
    if (allow_esac_as_next)
      {
	allow_esac_as_next--;
	if (STREQ (token, "esac"))
	  {
	    in_case_pattern_list = 0;
	    return (ESAC);
	  }
      }

    /* Ditto for `{' in the FUNCTION case. */
    if (allow_open_brace)
      {
	allow_open_brace = 0;
	if (token[0] == '{' && !token[1])
	  {
	    open_brace_awaiting_satisfaction++;
	    return ('{');
	  }
      }

    if (posixly_correct)
      CHECK_FOR_RESERVED_WORD (token);

#if defined (ALIAS)
    /* OK, we have a token.  Let's try to alias expand it, if (and only if)
       it's eligible. 

       It is eligible for expansion if the shell is in interactive mode, and
       the token is unquoted and the last token read was a command
       separator (or expand_next_token is set), and we are currently
       processing an alias (pushed_string_list is non-empty) and this
       token is not the same as the current or any previously
       processed alias.

       Special cases that disqualify:
	 In a pattern list in a case statement (in_case_pattern_list). */
    if (interactive_shell && !quoted && !in_case_pattern_list &&
	(expand_next_token || command_token_position (last_read_token)))
      {
	char *alias_expand_word (), *expanded;

	if (expanded_token_stack && token_has_been_expanded (token))
	  goto no_expansion;

	expanded = alias_expand_word (token);
	if (expanded)
	  {
	    int len = strlen (expanded), expand_next;

	    /* Erase the current token. */
	    token_index = 0;

	    expand_next = (expanded[len - 1] == ' ') ||
			  (expanded[len - 1] == '\t');

	    push_string (expanded, expand_next, token);
	    goto re_read_token;
	  }
	else
	  /* This is an eligible token that does not have an expansion. */
no_expansion:
	  expand_next_token = 0;
      }
    else
      {
	expand_next_token = 0;
      }
#endif /* ALIAS */

    if (!posixly_correct)
      CHECK_FOR_RESERVED_WORD (token);

    /* What if we are attempting to satisfy an open-brace grouper? */
    if (open_brace_awaiting_satisfaction && token[0] == '}' && !token[1])
      {
	open_brace_awaiting_satisfaction--;
	return ('}');
      }

    the_word = (WORD_DESC *)xmalloc (sizeof (WORD_DESC));
    the_word->word = xmalloc (1 + token_index);
    strcpy (the_word->word, token);
    the_word->dollar_present = dollar_present;
    the_word->quoted = quoted;
    the_word->assignment = assignment (token);

    yylval.word = the_word;
    result = WORD;

    /* A word is an assignment if it appears at the beginning of a
       simple command, or after another assignment word.  This is
       context-dependent, so it cannot be handled in the grammar. */
    if (assignment_acceptable (last_read_token) && the_word->assignment)
      result = ASSIGNMENT_WORD;

    if (last_read_token == FUNCTION)
      allow_open_brace = 1;
  }
  return (result);
}

/* Return 1 if TOKEN is a token that after being read would allow
   a reserved word to be seen, else 0. */
static int
reserved_word_acceptable (token)
     int token;
{
#if 0
  if (member (token, "\n;()|&{") ||
#else
  if (token == '\n' || token == ';' || token == '(' || token == ')' ||
      token == '|' || token == '&' || token == '{' ||
#endif
      token == '}' ||			/* XXX */
      token == AND_AND ||
      token == BANG ||
      token == DO ||
      token == ELIF ||
      token == ELSE ||
      token == FI ||
      token == IF ||
      token == OR_OR ||
      token == SEMI_SEMI ||
      token == THEN ||
      token == UNTIL ||
      token == WHILE ||
      token == DONE ||		/* XXX these two are experimental */
      token == ESAC ||
      token == 0)
    return (1);
  else
    return (0);
}

/* Return the index of TOKEN in the alist of reserved words, or -1 if
   TOKEN is not a shell reserved word. */
int
find_reserved_word (token)
     char *token;
{
  int i;
  for (i = 0; word_token_alist[i].word != (char *)NULL; i++)
    if (STREQ (token, word_token_alist[i].word))
      return i;
  return -1;
}

#if defined (READLINE)
/* Called after each time readline is called.  This insures that whatever
   the new prompt string is gets propagated to readline's local prompt
   variable. */
static void
reset_readline_prompt ()
{
  if (prompt_string_pointer)
    {
      char *temp_prompt;

      temp_prompt = *prompt_string_pointer
			? decode_prompt_string (*prompt_string_pointer)
			: (char *)NULL;

      if (temp_prompt == 0)
	{
	  temp_prompt = xmalloc (1);
	  temp_prompt[0] = '\0';
	}

      FREE (current_readline_prompt);

      current_readline_prompt = temp_prompt;
    }
}
#endif /* READLINE */

#if defined (HISTORY)
/* A list of tokens which can be followed by newlines, but not by
   semi-colons.  When concatenating multiple lines of history, the
   newline separator for such tokens is replaced with a space. */
static int no_semi_successors[] = {
  '\n', '{', '(', ')', ';', '&', '|',
  CASE, DO, ELSE, IF, IN, SEMI_SEMI, THEN, UNTIL, WHILE, AND_AND, OR_OR,
  0
};

/* If we are not within a delimited expression, try to be smart
   about which separators can be semi-colons and which must be
   newlines. */
char *
history_delimiting_chars ()
{
  if (!delimiter_depth)
    {
      register int i;

      for (i = 0; no_semi_successors[i]; i++)
	{
	  if (token_before_that == no_semi_successors[i])
	    return (" ");
	}
      return ("; ");
    }
  else
    return ("\n");
}
#endif /* HISTORY */

/* Issue a prompt, or prepare to issue a prompt when the next character
   is read. */
static void
prompt_again ()
{
  char *temp_prompt;

  if (!interactive)	/* XXX */
    return;

  ps1_prompt = get_string_value ("PS1");
  ps2_prompt = get_string_value ("PS2");

  if (!prompt_string_pointer)
    prompt_string_pointer = &ps1_prompt;

  temp_prompt = (*prompt_string_pointer)
			? decode_prompt_string (*prompt_string_pointer)
			: (char *)NULL;

  if (temp_prompt == 0)
    {
      temp_prompt = xmalloc (1);
      temp_prompt[0] = '\0';
    }

  current_prompt_string = *prompt_string_pointer;
  prompt_string_pointer = &ps2_prompt;

#if defined (READLINE)
  if (!no_line_editing)
    {
      FREE (current_readline_prompt);
      current_readline_prompt = temp_prompt;
    }
  else
#endif	/* READLINE */
    {
      FREE (current_decoded_prompt);
      current_decoded_prompt = temp_prompt;
    }
}

static void
print_prompt ()
{
  fprintf (stderr, "%s", current_decoded_prompt);
  fflush (stderr);
}

/* Return a string which will be printed as a prompt.  The string
   may contain special characters which are decoded as follows:
   
	\t	the time
	\d	the date
	\n	CRLF
	\s	the name of the shell
	\w	the current working directory
	\W	the last element of PWD
	\u	your username
	\h	the hostname
	\#	the command number of this command
	\!	the history number of this command
	\$	a $ or a # if you are root
	\<octal> character code in octal
	\\	a backslash
*/
#define PROMPT_GROWTH 50
char *
decode_prompt_string (string)
     char *string;
{
  int result_size = PROMPT_GROWTH;
  int result_index = 0;
  char *result;
  int c;
  char *temp = (char *)NULL;
  WORD_LIST *list;

#if defined (PROMPT_STRING_DECODE)

  result = xmalloc (PROMPT_GROWTH);
  result[0] = 0;

  while (c = *string++)
    {
      if (posixly_correct && c == '!')
	{
	  if (*string == '!')
	    {
	      temp = savestring ("!");
	      goto add_string;
	    }
	  else
	    {
#if !defined (HISTORY)
		temp = savestring ("1");
#else /* HISTORY */
		temp = itos (history_number ());
#endif /* HISTORY */
		string--;	/* add_string increments string again. */
		goto add_string;
	    }
	} 
      if (c == '\\')
	{
	  c = *string;

	  switch (c)
	    {
	    case '0':
	    case '1':
	    case '2':
	    case '3':
	    case '4':
	    case '5':
	    case '6':
	    case '7':
	      {
		char octal_string[4];
		int n;

		strncpy (octal_string, string, 3);
		octal_string[3] = '\0';

		n = read_octal (octal_string);
		temp = xmalloc (3);

		if (n == CTLESC || n == CTLNUL)
		  {
		    string += 3;
		    temp[0] = CTLESC;
		    temp[1] = n;
		    temp[2] = '\0';
		  }
		else if (n == -1)
		  {
		    temp[0] = '\\';
		    temp[1] = '\0';
		  }
		else
		  {
		    string += 3;
		    temp[0] = n;
		    temp[1] = '\0';
		  }

		c = 0;
		goto add_string;
	      }
	  
	    case 't':
	    case 'd':
	      /* Make the current time/date into a string. */
	      {
		time_t the_time = time (0);
		char *ttemp = ctime (&the_time);
		temp = savestring (ttemp);

		if (c == 't')
		  {
		    strcpy (temp, temp + 11);
		    temp[8] = '\0';
		  }
		else
		  temp[10] = '\0';

		goto add_string;
	      }

	    case 'n':
	      if (!no_line_editing)
		temp = savestring ("\r\n");
	      else
		temp = savestring ("\n");
	      goto add_string;

	    case 's':
	      {
		temp = base_pathname (shell_name);
		temp = savestring (temp);
		goto add_string;
	      }
	
	    case 'w':
	    case 'W':
	      {
		/* Use the value of PWD because it is much more effecient. */
#define EFFICIENT
#if defined(EFFICIENT)
		char *polite_directory_format (), t_string[MAXPATHLEN];

		temp = get_string_value ("PWD");

		if (!temp)
		  getwd (t_string);
		else
		  strcpy (t_string, temp);
#else
		getwd (t_string);
#endif	/* EFFICIENT */

		if (c == 'W')
		  {
		    char *dir = (char *)strrchr (t_string, '/');
		    if (dir && dir != t_string)
		      strcpy (t_string, dir + 1);
		    temp = savestring (t_string);
		  }
		else
		  temp = savestring (polite_directory_format (t_string));
		goto add_string;
	      }
      
	    case 'u':
	      {
		temp = savestring (current_user.user_name);
		goto add_string;
	      }

	    case 'h':
	      {
		char *t_string;

		temp = savestring (current_host_name);
		if (t_string = (char *)strchr (temp, '.'))
		  *t_string = '\0';
		goto add_string;
	      }

	    case '#':
	      {
		temp = itos (current_command_number);
		goto add_string;
	      }

	    case '!':
	      {
#if !defined (HISTORY)
		temp = savestring ("1");
#else /* HISTORY */
		temp = itos (history_number ());
#endif /* HISTORY */
		goto add_string;
	      }

	    case '$':
	      temp = savestring (geteuid () == 0 ? "#" : "$");
	      goto add_string;

#if defined (READLINE)
	    case '[':
	    case ']':
	      temp = xmalloc(3);
	      temp[0] = '\001';
	      temp[1] = (c == '[') ? RL_PROMPT_START_IGNORE : RL_PROMPT_END_IGNORE;
	      temp[2] = '\0';
	      goto add_string;
#endif

	    case '\\':
	      temp = savestring ("\\");
	      goto add_string;

	    default:
	      temp = savestring ("\\ ");
	      temp[1] = c;

	    add_string:
	      if (c)
		string++;
	      result =
		sub_append_string (temp, result, &result_index, &result_size);
	      temp = (char *)NULL; /* Free ()'ed in sub_append_string (). */
	      result[result_index] = '\0';
	      break;
	    }
	}
      else
	{
	  while (3 + result_index > result_size)
	    result = xrealloc (result, result_size += PROMPT_GROWTH);

	  result[result_index++] = c;
	  result[result_index] = '\0';
	}
    }
#else /* !PROMPT_STRING_DECODE */
  result = savestring (string);
#endif /* !PROMPT_STRING_DECODE */

  /* Perform variable and parameter expansion and command substitution on
     the prompt string. */
  list = expand_string_unsplit (result, 1);
  free (result);
  result = string_list (list);
  dispose_words (list);

  return (result);
}

/* Report a syntax error, and restart the parser.  Call here for fatal
   errors. */
yyerror ()
{
  report_syntax_error ((char *)NULL);
  reset_parser ();
}

/* Report a syntax error with line numbers, etc.
   Call here for recoverable errors.  If you have a message to print,
   then place it in MESSAGE, otherwise pass NULL and this will figure
   out an appropriate message for you. */
static void
report_syntax_error (message)
     char *message;
{
  if (message)
    {
      if (!interactive)
	{
	  char *name = bash_input.name ? bash_input.name : "stdin";
	  report_error ("%s: line %d: `%s'", name, line_number, message);
	}
      else
	{
	  if (EOF_Reached)
	    EOF_Reached = 0;
	  report_error ("%s", message);
	}

      last_command_exit_value = EX_USAGE;
      return;
    }

  if (shell_input_line && *shell_input_line)
    {
      char *t = shell_input_line;
      register int i = shell_input_line_index;
      int token_end = 0;

      if (!t[i] && i)
	i--;

      while (i && (t[i] == ' ' || t[i] == '\t' || t[i] == '\n'))
	i--;

      if (i)
	token_end = i + 1;

      while (i && !member (t[i], " \n\t;|&"))
	i--;

      while (i != token_end && member (t[i], " \t\n"))
	i++;

      if (token_end)
	{
	  char *error_token;
	  error_token = xmalloc (1 + (token_end - i));
	  strncpy (error_token, t + i, token_end - i);
	  error_token[token_end - i] = '\0';

	  report_error ("syntax error near unexpected token `%s'", error_token);
	  free (error_token);
	}
      else if ((i == 0) && (token_end == 0))	/* a 1-character token */
	{
	  char etoken[2];
	  etoken[0] = t[i];
	  etoken[1] = '\0';

	  report_error ("syntax error near unexpected token `%s'", etoken);
	}

      if (!interactive)
	{
	  char *temp = savestring (shell_input_line);
	  char *name = bash_input.name ? bash_input.name : "stdin";
	  int l = strlen (temp);

	  while (l && temp[l - 1] == '\n')
	    temp[--l] = '\0';

	  report_error ("%s: line %d: `%s'", name, line_number, temp);
	  free (temp);
	}
    }
  else
    {
      char *name, *msg;
      if (!interactive)
	name = bash_input.name ? bash_input.name : "stdin";
      if (EOF_Reached)
	msg = "syntax error: unexpected end of file";
      else
	msg = "syntax error";
      if (!interactive)
	report_error ("%s: line %d: %s", name, line_number, msg);
      else
	{
	  /* This file uses EOF_Reached only for error reporting
	     when the shell is interactive.  Other mechanisms are 
	     used to decide whether or not to exit. */
	  EOF_Reached = 0;
	  report_error (msg);
	}
    }
  last_command_exit_value = EX_USAGE;
}

/* ??? Needed function. ??? We have to be able to discard the constructs
   created during parsing.  In the case of error, we want to return
   allocated objects to the memory pool.  In the case of no error, we want
   to throw away the information about where the allocated objects live.
   (dispose_command () will actually free the command. */
discard_parser_constructs (error_p)
     int error_p;
{
}
   
/* Do that silly `type "bye" to exit' stuff.  You know, "ignoreeof". */

/* A flag denoting whether or not ignoreeof is set. */
int ignoreeof = 0;

/* The number of times that we have encountered an EOF character without
   another character intervening.  When this gets above the limit, the
   shell terminates. */
int eof_encountered = 0;

/* The limit for eof_encountered. */
int eof_encountered_limit = 10;

/* If we have EOF as the only input unit, this user wants to leave
   the shell.  If the shell is not interactive, then just leave.
   Otherwise, if ignoreeof is set, and we haven't done this the
   required number of times in a row, print a message. */
static void
handle_eof_input_unit ()
{
  if (interactive)
    {
      /* shell.c may use this to decide whether or not to write out the
	 history, among other things.  We use it only for error reporting
	 in this file. */
      if (EOF_Reached)
	EOF_Reached = 0;

      /* If the user wants to "ignore" eof, then let her do so, kind of. */
      if (ignoreeof)
	{
	  if (eof_encountered < eof_encountered_limit)
	    {
	      fprintf (stderr, "Use \"%s\" to leave the shell.\n",
		       login_shell ? "logout" : "exit");
	      eof_encountered++;
	      /* Reset the prompt string to be $PS1. */
	      prompt_string_pointer = (char **)NULL;
	      prompt_again ();
	      last_read_token = current_token = '\n';
	      return;
	    } 
	}

      /* In this case EOF should exit the shell.  Do it now. */
      reset_parser ();
      exit_builtin ((WORD_LIST *)NULL);
    }
  else
    {
      /* We don't write history files, etc., for non-interactive shells. */
      EOF_Reached = 1;
    }
}
