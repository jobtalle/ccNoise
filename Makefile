# Generic makefile for static libraries

NAME=ccNoise

SOURCEDIR=src/$(NAME)
LIBDIR=lib
INCDIR=include
BINDIR=bin
TESTDIR=test

RM=rm -f
AR=ar rcs
CFLAGS=-I$(INCDIR) -O3 -DCC_USE_ALL
LDLIBS=-lGL -lGLU -lGLEW -lm

SRCS=$(wildcard ./$(SOURCEDIR)/*.c)
OBJS=$(subst .c,.o,$(SRCS))
LIBFILE=lib$(NAME).a
MAKEFILEDIR=$(dir $(abspath $(lastword $(MAKEFILE_LIST))))

all: $(NAME)

.PHONY: $(NAME)
$(NAME): $(OBJS)
	$(AR) $(LIBDIR)/lib$(NAME).a $(OBJS)

.PHONY: test
test: $(NAME)
	@(cd $(TESTDIR); $(MAKE) BINDIR="$(MAKEFILEDIR)$(BINDIR)" INCDIR="$(MAKEFILEDIR)$(INCDIR)" LIBDIR="$(MAKEFILEDIR)$(LIBDIR)" LIBNAME="$(NAME)")

.PHONY: clean
clean:
	find $(SOURCEDIR) -type f -name '*.o' -delete
	$(RM) $(LIBDIR)/$(LIBFILE)

.PHONY: install
install:
	mkdir -p $(DESTDIR)/usr/include
	cp -R $(INCDIR)/* $(DESTDIR)/usr/include
	mkdir -p $(DESTDIR)/usr/lib
	cp -R $(LIBDIR)/* $(DESTDIR)/usr/lib
	@(cd $(UTILDIR); $(MAKE) install BINDIR="$(MAKEFILEDIR)$(BINDIR)" INCDIR="$(MAKEFILEDIR)$(INCDIR)" LIBDIR="$(MAKEFILEDIR)$(LIBDIR)" LIBNAME="$(NAME)")

.PHONY: dist-clean
dist-clean: clean
	$(RM) *~ .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CC) $(CFLAGS) -MM $(SRCS) >>./.depend;

include .depend
