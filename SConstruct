import os

env=Environment(CC='gcc')

sources=[Glob('test/*.c')]
libs=['ccore', 'X11', 'Xrandr', 'Xinerama', 'Xi', 'GL', 'GLU', 'GLEW', 'pthread', 'ccNoise', 'ccTrigonometry', 'ccRandom', 'ccSort', 'm']

env.Append(CCFLAGS=['-DCC_USE_ALL'])

env.Append(CCFLAGS=['-g'])
env.Append(CCFLAGS=['-Wall'])

env.Program(target='bin/test', source=sources, LIBS=libs)
