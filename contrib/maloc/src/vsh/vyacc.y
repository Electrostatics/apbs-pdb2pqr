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
 * rcsid="$Id: vyacc.y,v 1.7 2008/03/12 05:13:59 fetk Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     vyacc.y
 *
 * Purpose:  Class Vsh: YACC Specification
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

%{

#include "vsh_p.h"

VEMBED(rcsid="$Id: vyacc.y,v 1.7 2008/03/12 05:13:59 fetk Exp $")

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

%}

%union {
    WORD_DESC *word;         /* the word that we read. */
    /* int number; */        /* number we saw */
    WORD_LIST *word_list;    /* a sequence of white-space separated words */
    COMMAND *command;        /* a complete command */
    REDIRECT *redirect;      /* redirect i/o info */
    ELEMENT element;         /* base element */
    PATTERN_LIST *pattern;   /* a case pattern */
}

/* Reserved words coming from yylex() */
%token IF THEN ELSE ELIF FI
%token CASE IN ESAC
%token FOR WHILE UNTIL DO DONE

/* Special 2-character tokens coming from yylex() */
%token SEMI_SEMI AND_AND OR_OR LESS_LESS GREATER_GREATER
%token LESS_AND GREATER_AND AND_GREATER LESS_GREATER GREATER_BAR

/* Special 3-character tokens coming from yylex() */
%token LESS_LESS_MINUS

/* Remaining tokens coming from yylex() */
%token <word> WORD ASSIGNMENT_WORD

%token IF THEN ELSE ELIF FI CASE ESAC FOR SELECT WHILE UNTIL DO DONE FUNCTION
%token IN BANG

/* New Non-terminal types */
%type <command> for_command case_command while_command until_command

/* Non-terminal types */
%type <command> inputunit command pipeline
%type <command> list list0 list1 simple_list simple_list1
%type <command> simple_command shell_command_1 shell_command
%type <command> group_command function_def if_command elif_clause subshell
%type <redirect> redirection redirections
%type <element> simple_command_element
%type <word_list> words pattern
%type <pattern> pattern_list case_clause_sequence case_clause_1 pattern_list_1

/* Start rule */
%start inputunit

/* Precedence (low-to-high) */
%left '&' ';' '\n' vsh_EOF
%left AND_AND OR_OR
%right '|'

/* The complete grammar specification */
%%

   /*
    * ************************************************************************
    * start state
    * ************************************************************************
    */

inputunit : simple_list '\n' {
                Vsh_trace("Yacc","inputunit::1");
                global_command = $1;
                if (global_command != VNULL) {
                    if (global_command->type == cm_simple) {
                        wordTmp = $1->value.Simple->words;
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
          | '\n' {
                Vsh_trace("Yacc","inputunit::2");
                global_command = (COMMAND*)VNULL;
                YYACCEPT;
            }
          | vsh_EOF {
                Vsh_trace("Yacc","inputunit::3");
                global_command = (COMMAND*)VNULL;
                cmdKey = 2;
                YYACCEPT;
            }
          | error newlines {
                Vsh_trace("Yacc","inputunit::4");
                global_command = (COMMAND*)VNULL;
                YYABORT;
            }
          ;

   /*
    * ************************************************************************
    * sim lists and regular lists
    * ************************************************************************
    */

   /*
    * A simple_list is a list that contains no significant newline
    * and no leading or trailing newline.  newline are allowed
    * only following operators, where they are not significant.
    * This is what an inputunit consists of.
    */
simple_list : simple_list1 {
                  Vsh_trace("Yacc","simple_list::1");
                  $$ = $1;
              }
            | simple_list1 '&' {
                  Vsh_trace("Yacc","simple_list::2");
                  $$ = $1;
              }
            | simple_list1 ';' {
                  Vsh_trace("Yacc","simple_list::3");
                  $$ = $1;
              }
            ;

simple_list1 : simple_list1 AND_AND newlines simple_list1 {
                   Vsh_trace("Yacc","simple_list1::1");
                   $$ = $1;
               }
             | simple_list1 OR_OR newlines simple_list1 {
                   Vsh_trace("Yacc","simple_list1::2");
                   $$ = $1;
               }
             | simple_list1 '&' simple_list1 {
                   Vsh_trace("Yacc","simple_list1::3");
                   $$ = $1;
               }
             | simple_list1 ';' simple_list1 {
                   Vsh_trace("Yacc","simple_list1::4");
                   $$ = $1;
               }
             | pipeline {
                   Vsh_trace("Yacc","simple_list1::5");
                   $$ = $1;
               }
             | '!' pipeline {
                   Vsh_trace("Yacc","simple_list1::6");
                   $$ = $2;
               }
             ;

pipeline : pipeline '|' newlines pipeline {
               Vsh_trace("Yacc","pipeline::1");
               $$ = $1;
           }
         | command {
               Vsh_trace("Yacc","pipeline::2");
               $$ = $1;
           }
         ;

   /*
    * A list allows leading or trailing newline and
    * newline as operators (equivalent to semicolons).
    * It must end with a newline or semicolon.
    * Lists are used within commands such as if, for, while.
    */
list : newlines list0 {
           Vsh_trace("Yacc","list::1");
           $$ = $2;
       }
     ;

list0 : list1 {
            Vsh_trace("Yacc","list0::1");
            $$ = $1;
        }
      | list1 '\n' newlines {
            Vsh_trace("Yacc","list0::2");
            $$ = $1;
        }
      | list1 '&' newlines {
            Vsh_trace("Yacc","list0::3");
            $$ = $1;
        }
      | list1 ';' newlines {
            Vsh_trace("Yacc","list0::4");
            $$ = $1;
        }
      ;

list1 : list1 AND_AND newlines list1 {
            Vsh_trace("Yacc","list1::1");
            $$ = $1;
        }
      | list1 OR_OR newlines list1 {
            Vsh_trace("Yacc","list1::2");
            $$ = $1;
        }
      | list1 '&' newlines list1 {
            Vsh_trace("Yacc","list1::3");
            $$ = $1;
        }
      | list1 ';' newlines list1 {
            Vsh_trace("Yacc","list1::4");
            $$ = $1;
        }
      | list1 '\n' newlines list1 {
            Vsh_trace("Yacc","list1::5");
            $$ = $1;
        }
      | pipeline {
            Vsh_trace("Yacc","list1::6");
            $$ = $1;
        }
      | '!' pipeline {
            Vsh_trace("Yacc","list1::7");
            $$ = $2;
        }
      ;

   /*
    * ************************************************************************
    * redirection (creates MANY shift/reduce conflicts unfortunately...)
    * ************************************************************************
    */

redirection : '>' WORD {
                  Vsh_trace("Yacc","redirection::1");
                  redir.filename = $2;
                  $$ = make_redirection(1, r_output_direction, redir);
              }
            | GREATER_GREATER WORD {
                  Vsh_trace("Yacc","redirection::2");
                  redir.filename = $2;
                  $$ = make_redirection(1, r_appending_to, redir);
              }

           /*
            : '>' WORD { Vsh_trace("Yacc","redirection"); }
            | '<' WORD { Vsh_trace("Yacc","redirection"); }
            | GREATER_GREATER WORD { Vsh_trace("Yacc","redirection"); }
            | LESS_LESS WORD { Vsh_trace("Yacc","redirection"); }
            | number '>' WORD { Vsh_trace("Yacc","redirection"); }
            | number '<' WORD { Vsh_trace("Yacc","redirection"); }
            | number GREATER_GREATER WORD { Vsh_trace("Yacc","redirection"); }
            | number LESS_LESS WORD { Vsh_trace("Yacc","redirection"); }
            | LESS_AND number { Vsh_trace("Yacc","redirection"); }
            | number LESS_AND number { Vsh_trace("Yacc","redirection"); }
            | GREATER_AND number { Vsh_trace("Yacc","redirection"); }
            | number GREATER_AND number { Vsh_trace("Yacc","redirection"); }
            | number LESS_AND WORD { Vsh_trace("Yacc","redirection"); }
            | number GREATER_AND WORD { Vsh_trace("Yacc","redirection"); }
            | number LESS_LESS_MINUS WORD { Vsh_trace("Yacc","redirection"); }
            | number GREATER_AND '-' { Vsh_trace("Yacc","redirection"); }
            | number LESS_AND '-' { Vsh_trace("Yacc","redirection"); }
            | number LESS_GREATER WORD { Vsh_trace("Yacc","redirection"); }
            | number GREATER_BAR WORD { Vsh_trace("Yacc","redirection"); }
            | GREATER_GREATER WORD { Vsh_trace("Yacc","redirection"); }
            | LESS_LESS WORD { Vsh_trace("Yacc","redirection"); }
            | LESS_AND WORD { Vsh_trace("Yacc","redirection"); }
            | GREATER_AND WORD { Vsh_trace("Yacc","redirection"); }
            | LESS_LESS_MINUS WORD { Vsh_trace("Yacc","redirection"); }
            | GREATER_AND '-' { Vsh_trace("Yacc","redirection"); }
            | LESS_AND '-' { Vsh_trace("Yacc","redirection"); }
            | AND_GREATER WORD { Vsh_trace("Yacc","redirection"); }
            | LESS_GREATER WORD { Vsh_trace("Yacc","redirection"); }              
            | GREATER_BAR WORD { Vsh_trace("Yacc","redirection"); }
           */

            ;

redirections : redirection {
                   Vsh_trace("Yacc","redirections::1");
                   $$ = $1;
               }
             | redirections redirection {
                   Vsh_trace("Yacc","redirections::2");
                   $$ = $1;
               }
             ;

   /*
    * ************************************************************************
    * sim and shell commands
    * ************************************************************************
    */

command : simple_command {
              Vsh_trace("Yacc","command::1");
              $$ = clean_simple_command($1);
          }
        | shell_command {
              Vsh_trace("Yacc","command::2");
              $$ = $1;
          }
        ;

simple_command : simple_command_element {
                     Vsh_trace("Yacc","simple_command::1");
                     $$ = make_simple_command($1,(COMMAND*)VNULL);
                 }
               | simple_command simple_command_element {
                     Vsh_trace("Yacc","simple_command::2");
                     $$ = make_simple_command($2,$1);
                 }
               ;

shell_command : shell_command_1 {
                    Vsh_trace("Yacc","shell_command::1");
                    $$ = $1;
                }
              | shell_command_1 redirections {
                    Vsh_trace("Yacc","shell_command::2");
                    $$ = $1;
                }
              ;

shell_command_1 : for_command {
                      Vsh_trace("Yacc","shell_command_1::1");
                      $$ = VNULL;
                  }
                | case_command {
                      Vsh_trace("Yacc","shell_command_1::2");
                      $$ = VNULL;
                  }
                | while_command {
                      Vsh_trace("Yacc","shell_command_1::3");
                      $$ = VNULL;
                  }
                | until_command {
                      Vsh_trace("Yacc","shell_command_1::4");
                      $$ = VNULL;
                  }
                | if_command {
                      Vsh_trace("Yacc","shell_command_1::5");
                      $$ = VNULL;
                  }
                | subshell {
                      Vsh_trace("Yacc","shell_command_1::6");
                      $$ = VNULL;
                  }
                | group_command {
                      Vsh_trace("Yacc","shell_command_1::7");
                      $$ = VNULL;
                  }
                | function_def {
                      Vsh_trace("Yacc","shell_command_1::8");
                      $$ = VNULL;
                  }
                ;

simple_command_element : WORD {
                             Vsh_trace("Yacc","simple_command_element::1");
                             $$.word     = $1;
                             $$.redirect = 0;
                         }
                       | ASSIGNMENT_WORD {
                             Vsh_trace("Yacc","simple_command_element::2");
                             $$.word     = $1;
                             $$.redirect = 0;
                         }
                       | redirection {
                             Vsh_trace("Yacc","simple_command_element::3");
                             $$.word     = 0;
                             $$.redirect = $1;
                         }
                       ;

   /*
    * ************************************************************************
    * subshell, group, function definition
    * ************************************************************************
    */

subshell : '(' list ')' {
               Vsh_trace("Yacc","subshell::1");
               $$ = $2;
           }
         ;
    
group_command : '{' list '}' {
                    Vsh_trace("Yacc","group_command::1");
                    $$ = $2;
                }
              ;

function_def : WORD '(' ')' newlines group_command {
                   Vsh_trace("Yacc","function_def::1");
                   $$ = VNULL;
               }
             | WORD '(' ')' newlines group_command redirections {
                   Vsh_trace("Yacc","function_def::2");
                   $$ = VNULL;
               }
             ;

   /*
    * ************************************************************************
    * if command
    * ************************************************************************
    */

if_command : IF list THEN list FI {
                 Vsh_trace("Yacc","if_command::1");
                 $$ = VNULL;
             }
           | IF list THEN list ELSE list FI {
                 Vsh_trace("Yacc","if_command::2");
                 $$ = VNULL;
             }
           | IF list THEN list elif_clause FI {
                 Vsh_trace("Yacc","if_command::3");
                 $$ = VNULL;
             }
           ;

elif_clause : ELIF list THEN list {
                  Vsh_trace("Yacc","elif_clause::1");
                  $$ = VNULL;
              }
            | ELIF list THEN list ELSE list {
                  Vsh_trace("Yacc","elif_clause::2");
                  $$ = VNULL;
              }
            | ELIF list THEN list elif_clause {
                  Vsh_trace("Yacc","elif_clause::3");
                  $$ = VNULL;
              }
            ;

   /*
    * ************************************************************************
    * case command
    * ************************************************************************
    */

case_command : CASE WORD newlines IN newlines ESAC {
                   Vsh_trace("Yacc","case_command::1");
                   $$ = VNULL;
               }
             | CASE WORD newlines IN case_clause_sequence newlines ESAC {
                   Vsh_trace("Yacc","case_command::2");
                   $$ = VNULL;
               }
             | CASE WORD newlines IN case_clause_1 ESAC {
                   Vsh_trace("Yacc","case_command::3");
                   $$ = VNULL;
               }
             ;

case_clause_1 : pattern_list_1 {
                    Vsh_trace("Yacc","case_clause_1::1");
                    $$ = VNULL;
                }
              | case_clause_sequence pattern_list_1 {
                    Vsh_trace("Yacc","case_clause_1::2");
                    $$ = VNULL;
                }
              ;

pattern_list_1 : newlines pattern ')' list {
                     Vsh_trace("Yacc","pattern_list_1::1");
                     $$ = VNULL;
                 }
               | newlines pattern ')' newlines {
                     Vsh_trace("Yacc","pattern_list_1::2");
                     $$ = VNULL;
                 }
               | newlines '(' pattern ')' list {
                     Vsh_trace("Yacc","pattern_list_1::3");
                     $$ = VNULL;
                 }
               | newlines '(' pattern ')' newlines {
                     Vsh_trace("Yacc","pattern_list_1::4");
                     $$ = VNULL;
                 }
               ;

case_clause_sequence : pattern_list {
                           Vsh_trace("Yacc","case_clause_sequence::1");
                           $$ = VNULL;
                       }
                     | case_clause_sequence pattern_list {
                           Vsh_trace("Yacc","case_clause_sequence::2");
                           $$ = VNULL;
                       }
                     ;

pattern_list : newlines pattern ')' list SEMI_SEMI {
                   Vsh_trace("Yacc","pattern_list::1");
                   $$ = VNULL;
               }
             | newlines pattern ')' newlines SEMI_SEMI {
                   Vsh_trace("Yacc","pattern_list::2");
                   $$ = VNULL;
               }
             | newlines '(' pattern ')' list SEMI_SEMI {
                   Vsh_trace("Yacc","pattern_list::3");
                   $$ = VNULL;
               }
             | newlines '(' pattern ')' newlines SEMI_SEMI {
                   Vsh_trace("Yacc","pattern_list::4");
                   $$ = VNULL;
               }
             ;

pattern : WORD {
              Vsh_trace("Yacc","pattern::1");
              $$ = VNULL;
          }
        | pattern '|' WORD {
              Vsh_trace("Yacc","pattern::2");
              $$ = VNULL;
          }
        ;

   /*
    * ************************************************************************
    * for command
    * ************************************************************************
    */

for_command : FOR WORD newlines DO list DONE {
                  Vsh_trace("Yacc","for_command::1");
                  $$ = VNULL;
              }
            | FOR WORD newlines '{' list '}' {
                  Vsh_trace("Yacc","for_command::2");
                  $$ = VNULL;
              }
            | FOR WORD ';' newlines DO list DONE {
                  Vsh_trace("Yacc","for_command::3");
                  $$ = VNULL;
              }
            | FOR WORD ';' newlines '{' list '}' {
                  Vsh_trace("Yacc","for_command::4");
                  $$ = VNULL;
              }
            | FOR WORD newlines IN words term newlines DO list DONE {
                  Vsh_trace("Yacc","for_command::5");
                  $$ = VNULL;
              }
            | FOR WORD newlines IN words term newlines '{' list '}' {
                  Vsh_trace("Yacc","for_command::6");
                  $$ = VNULL;
              }
            ;

   /*
    * ************************************************************************
    * while command
    * ************************************************************************
    */

while_command : WHILE list DO list DONE {
                    Vsh_trace("Yacc","while_command::1");
                    $$ = VNULL;
                }
              ;

   /*
    * ************************************************************************
    * until command
    * ************************************************************************
    */

until_command : UNTIL list DO list DONE {
                    Vsh_trace("Yacc","until_command::1");
                    $$ = VNULL;
                }
              ;

   /*
    * ************************************************************************
    * word, words, term, newline
    * ************************************************************************
    */

words : {
            Vsh_trace("Yacc","words::1");
            $$ = (WORD_LIST*)VNULL;
        }
      | words WORD {
            Vsh_trace("Yacc","words::2");
            $$ = make_word_list($2,$1);
        }
      ;

term : '\n' {
           Vsh_trace("Yacc","term::1");
       }
     | ';' {
           Vsh_trace("Yacc","term::2");
       }
     | vsh_EOF {
           Vsh_trace("Yacc","term::3");
       }
     ;

newlines : {
              Vsh_trace("Yacc","newlines::1");
          }
        | newlines '\n' {
              Vsh_trace("Yacc","newlines::2");
          }
        ;

%%

