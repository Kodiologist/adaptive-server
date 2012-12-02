adaptive-server is a bunch of R code that you can run under `Rserve`_ to provide a server that can fit models of intertemporal choice adaptively. Each time a human subject makes a choice, the choice is sent to the server, the server fits a Bayesian model of all the subject's choices so far (using either MCMC via `Stan`_ or grid approximation, depending on the model), and the server sends back a new question to ask the subject that in some sense maximizes the question's diagnostic value for model-fitting. You could run adaptive-server on the same machine as the client program, but the point of writing it as a server was that I could run it on a machine with lots of CPU and RAM while keeping my web server on a tiny VPS.

I wrote adaptive-server for use with `Builder`_.

Usage
============================================================

To start the server, try the shell command::

    R CMD Rserve --RS-conf Rserv.conf

``Rserv.conf`` should be a file containing something like::

    encoding utf8
    eval db.path = "adaptive-server-jobs.sqlite"
    source init-serv.R

The SQLite database should be initialized with the included ``jobs-schema.sql``.

License
============================================================

This program is copyright 2012 Kodi Arfer.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the `GNU General Public License`_ for more details.

.. _Rserve: http://www.rforge.net/Rserve/
.. _Stan: http://mc-stan.org
.. _Builder: https://github.com/Kodiologist/Builder
.. _`GNU General Public License`: http://www.gnu.org/licenses/
