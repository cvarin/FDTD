/* AUTORIGHTS
Copyright (C) 2006-2012  Ilker R. Capoglu

    This file is part of the Angora package.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//functions for handling file path strings

#include "headers.h"

#include "globals.h"

#include <cstring>
#include <cerrno>

//For file-directory manipulations
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <libgen.h>
#endif


void add_slash_to_path(string& path)
{//checks if there is a slash at the end of the path, and appends it if there is none
	if (path!="")
	{
		if (path[path.size()-1]!='/')
		{
			path = path + "/";
		}
	}
}

bool is_absolute_path(string& path)
{//checks if the path is absolute, by checking for the "/" character up front
	if (path!="")
	{
		return (path[0]=='/');
	}
	else
	{
		return false;
	}
}

#ifndef _WIN32
/* Function with behaviour like `mkdir -p'  */
// from the website "Niall's Weblog: The website of Niall Higgins (http://niallohiggins.com/2009/01/08/mkpath-mkdir-p-alike-in-c-for-unix/)"
// This function is specific to GNU
int mkpath(const char *s, mode_t mode){
        char *q, *r = NULL, *path = NULL, *up = NULL;
        int rv;

        rv = -1;
        if (strcmp(s, ".") == 0 || strcmp(s, "/") == 0)
                return (0);

        if ((path = strdup(s)) == NULL)
                exit(1);

        if ((q = strdup(s)) == NULL)
                exit(1);

        if ((r = dirname(q)) == NULL)
                goto out;

        if ((up = strdup(r)) == NULL)
                exit(1);

        if ((mkpath(up, mode) == -1) && (errno != EEXIST))
                goto out;

        if ((mkdir(path, mode) == -1) && (errno != EEXIST))
                rv = -1;
        else
                rv = 0;

out:
        if (up != NULL)
                free(up);
        free(q);
        free(path);
        return (rv);
}
#endif

int create_path(const string& path)
{//creates the given path, only if check mode is not enabled
#ifdef _WIN32
	int rv = mkdir(path.c_str());
#else
	int rv = mkpath(path.c_str(),(S_IRWXU|S_IRGRP|S_IROTH));
#endif
	return rv;  //<0 if not successful
}
