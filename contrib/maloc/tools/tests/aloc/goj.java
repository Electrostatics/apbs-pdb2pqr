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
 * rcsid="$Id: goj.java,v 1.5 2008/03/12 05:14:00 fetk Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     goj.java
 *
 * Purpose:  Test Java interaction with the MALOC library.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

import java.net.*;
import java.io.*;

public class goj {

    public static void main(String[] args) {

        int VPORTNUMBER = 14916; /* portbase;  5000 < VPORTNUMBER < 49152 */

        Socket theSocket;
        String hostname;
        PrintStream theOutputStream;

        hostname = "localhost";
        try {
            System.out.print("Connecting socket..");
            theSocket = new Socket(hostname, VPORTNUMBER + 1);
            theOutputStream = new PrintStream(theSocket.getOutputStream());
            System.out.print("..sending data..");

            theOutputStream.println("bhsingle");
            theOutputStream.println("putl 3");
            theOutputStream.println("0.0 1.0 0.0");
            theOutputStream.println("0.0 0.0 1.0");
            theOutputStream.println("0.0 0.0 0.0");
            theOutputStream.println("list 7");
            theOutputStream.println("fill 2 3");
            theOutputStream.println("0.0 1.0 0.0");
            theOutputStream.println("0.0 0.0 1.0");
            theOutputStream.println("0.5 0.5 0.5");
            theOutputStream.println("fill 3 3");
            theOutputStream.println("1.0 1.0 0.0");
            theOutputStream.println("0.0 1.0 1.0");
            theOutputStream.println("0.5 0.5 0.5");
            theOutputStream.println("list -7");
            theOutputStream.println("list 5");
            theOutputStream.println("line 1 4");
            theOutputStream.println("0.0 1.0 0.0");
            theOutputStream.println("0.0 0.0 1.0");
            theOutputStream.println("0.5 0.5 0.5");
            theOutputStream.println("0.0 1.0 0.0");
            theOutputStream.println("line 1 4");
            theOutputStream.println("1.0 1.0 0.0");
            theOutputStream.println("0.0 1.0 1.0");
            theOutputStream.println("0.5 0.5 0.5");
            theOutputStream.println("1.0 1.0 0.0");
            theOutputStream.println("list -5");
            theOutputStream.println("list 2");
            theOutputStream.println("fill 2 3");
            theOutputStream.println("0.0 1.0 0.0");
            theOutputStream.println("0.0 0.0 1.0");
            theOutputStream.println("0.5 0.5 0.5");
            theOutputStream.println("list -2");
            theOutputStream.println("list 3");
            theOutputStream.println("fill 3 3");
            theOutputStream.println("1.0 1.0 0.0");
            theOutputStream.println("0.0 1.0 1.0");
            theOutputStream.println("0.5 0.5 0.5");
            theOutputStream.println("list -3");
            theOutputStream.println("putl -1");

            System.out.print("..closing socket..");
            theSocket.close();
            System.out.print("..done.\n");
        }
        catch (UnknownHostException e) {
            System.err.println(e);
        }
        catch (IOException e) {
            System.err.println(e);
        }
    }
}

