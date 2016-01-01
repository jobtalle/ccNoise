import os

# Build library
env=Environment(CC='gcc', CCFLAGS='-Iinclude/')

sources=[Glob('src/*/*.c')]
libs=['GL', 'GLU', 'GLEW']

staticLibrary=env.Library(target='lib/ccNoise', source=sources, LIBS=libs)

# Build test
sources=[Glob('test/*.c')]
libs=['ccore', 'X11', 'Xrandr', 'Xinerama', 'Xi', 'GL', 'GLU', 'GLEW', 'pthread', 'm']
libpaths=['/usr/lib', '/usr/local/lib', '.']

env.Append(CCFLAGS=['-DCC_USE_ALL'])

env.Append(CCFLAGS=['-g'])
env.Append(CCFLAGS=['-Wall'])

env.Program(target='bin/test', source=sources, LIBS=[staticLibrary, libs], LIBPATH=libpaths)
