NAME=test

SOURCEDIR=.

RM=rm -f
CFLAGS:=-I$(INCDIR) -g -Wall -D_DEBUG -DCC_USE_ALL 
LDLIBS:=-L$(LIBDIR) $(foreach lib,$(LIBNAMES),-l$(lib)) -lccore -lGL -lGLU -lGLEW -lm -lX11 -lXrandr -lXinerama -lXi -lpthread

SRCS:=$(shell find $(SOURCEDIR) -name '*.c')
OBJS:=$(subst .c,.o,$(SRCS))

all: $(NAME)

.PHONY: $(NAME)
$(NAME): $(OBJS) .depend
	@(mkdir -p $(BINDIR))
	$(CC) $(LDFLAGS) -o $(BINDIR)/$(NAME) $(OBJS) $(LDLIBS)
	#@($(BINDIR)/$(NAME))

.PHONY: clean
clean:
	$(RM) $(OBJS)

.PHONY: dist-clean
dist-clean: clean
	$(RM) *~ .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CC) $(CFLAGS) -MM $(SRCS) >>./.depend;

include .depend
