pl034-1
=======

Scripts updated for 0.63.1.

loglm is now operational (needs loglin from 0.63.2 pre-release to work).

Minor corrections to addterm/droptem, lda and qda.

All the examples on man pages should now work in R, and there are now
many more of them.


pl030-1
=======

Scripts updated for 0.63.

stepAIC, addterm, dropterm now make use of add1 and drop1,  so
require 0.63.

More changes to use TRUE and FALSE instead of T and F, including in
help pages.

Checking for WIN32 in C code is replaced by checking for S-PLUS 4.x, so
the code will compile on Guido Masarotto's win32 port.


pl027-1
=======

plot.profile, pairs.profile and profile.glm are now implemented.

Column name in michelson is corrected to `Run'.

TRUE and FALSE used instead of T and F.

Additional functions, including addterm, dropterm.
