# Generic makefile for static libraries

NAME=ccNoise

SOURCEDIR=src
LIBDIR=lib
INCDIR=include
BINDIR=bin
TESTDIR=test

RM=rm -f
AR=ar rcs
CFLAGS=-I$(INCDIR) -O3 -DCC_USE_ALL
LDLIBS=-lGL -lGLU -lGLEW -lm

LIBS:=$(subst ./$(SOURCEDIR)/,,$(wildcard ./$(SOURCEDIR)/*))
SRCS:=$(wildcard ./$(SOURCEDIR)/*/*.c)
OBJS:=$(subst .c,.o,$(SRCS))
LIBFILE:=lib$(NAME).a
MAKEFILEDIR:=$(dir $(abspath $(lastword $(MAKEFILE_LIST))))

all: $(LIBS)

.PHONY: $(NAME)
$(LIBS): $(OBJS)
	$(AR) $(LIBDIR)/lib$@.a $(wildcard ./$(SOURCEDIR)/$@/*.o)

.PHONY: test
test: $(NAME)
	@(cd $(TESTDIR); $(MAKE) BINDIR="$(MAKEFILEDIR)$(BINDIR)" INCDIR="$(MAKEFILEDIR)$(INCDIR)" LIBDIR="$(MAKEFILEDIR)$(LIBDIR)" LIBNAMES="$(LIBS)")

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

.PHONY: dist-clean
dist-clean: clean
	$(RM) *~ .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CC) $(CFLAGS) -MM $(SRCS) >>./.depend;

include .depend
