#!/usr/bin/env python
# PyINDI: access INDI properties
# typical interpreter usage:
#   from pyindi import *
#   pi = PyINDI()
#
# Using onesession works, just beware overlapped usage such as in getFITS and that a
# long-pending open socket can get too full while idle.


import subprocess
from subprocess import Popen, PIPE
import sys
import numpy
import pyfits
import fnmatch
import socket
import time

class PyINDIException (Exception):
	"""Exception raised by PyINDI to report trouble
	"""
	pass

class PyINDI:
	"""Provide Python interaction with INDI properties.
	Requires getINDI, setINDI and evalINDI to be in execution PATH.
	Trouble is reported by raising exception PyINDIException.
	"""

	# local helper functions

	# Insure we have access to the given command
	def _trycmd (self, cmd):
	    try:
		# try to run the given command
		subprocess.check_call ([cmd, "--help"], stderr=PIPE) 
	    except OSError:
		# raised if command is not found at all
		raise PyINDIException("Can not find required program " + cmd)
	    except subprocess.CalledProcessError:
		# raised if command exits with other than status 0
		pass


	# Open a FITS file with the given BLOB name, with either fits or fts suffix
	def _openFITS (self, blob):
	    try:
		return (pyfits.open(blob+".fits"))
	    except IOError:
		try:
		    return (pyfits.open(blob+".fts"))
		except IOError:
		    raise PyINDIException(blob + " not found")

	# called when new instance created
	def __init__ (self, host="lbti-web", port=7624, verbose=False, onesession=False):
	    """Initialization for a new instance of PyINDI

	    Optional named arguments:
	        host: network host connection, default is lbti-web
	        port: network port connection, default is 7624
	        verbose:  whether to trace function execution, default False
		oneseesion: whether to open and reuse one connection for each command
	    """

	    # gather options
	    self.verbose = verbose
	    self.host = host
	    self.port = str(port)
	    self.onesession = onesession

	    # make sure we have access to get/set/evalINDI
	    self._trycmd ("getINDI")
	    self._trycmd ("setINDI")
	    self._trycmd ("evalINDI")

	    if self.onesession:
		# open the host/port for all subsequent commands
		self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
		self.s.connect((self.host, int(self.port)))
		self.fd = str(self.s.fileno())

	    if self.verbose:
	    	print "new PyINDI instance is ready"

	# fetch INDI properties
	def getINDI(self, *args, **kwargs):
	    """Fetch INDI properties

	    Optional named arguments:
		timeout: max time to wait for given properties, seconds; default 2
		showlabels: just display Properties with short descriptions,
		    nothing is returned; default False

	    Usage:
		Call with zero or more property name strings or lists of same:
		    return dictionary of each matching name:value
		Special case 1: if called with exactly one complete name:
		    return value
		Special case 2: if called with showlabels=True:
		    just display Properties with short descriptions, nothing returned

		Examples:
		    dict = getINDI ()			# returns all properties
		    dict = getINDI ("*.*.*")		# same with all wildcards
		    dict = getINDI (showlabels=True)	# display all properties with labels
		    dict = getINDI ("A.B.C", "D.E.F")	# call with two name strings
		    dict = getINDI (["A.B.C", "D.E.F"])	# call with one list of names
		    dict = getINDI ("A.B.*")		# use wildcard
		    valu = getINDI ("A.B.C")		# returns value, not dict
	    """

	    showlabels = kwargs.pop('showlabels', False)

	    # create a dictionary from a string of the form "A=1\nB=2\nC=3\n".
	    # allow for multiple lines in one value
	    def _crackGIP (s):
		d = {}
		for l in s.splitlines():
		    p = l.split('=')
		    if len(p) == 2:
			d[p[0]] = p[1]
			lastp = p[0]
		    else:
			d[lastp] += "\n" + p[0]
		return d

	    # convert INDI strig to native type
	    def _giType(s):
		if s == "On":
		    s = True
		elif s == "Off":
		    s = False
		else:
		    try:
			s = float("%.12g" % float(s))
		    except ValueError:
			pass
		return s

	    # initial command line
	    cmd = ["getINDI", "-h", self.host, "-p", self.port, "-w", "-t", str(kwargs.pop('timeout', 2))]
	    if showlabels:
		cmd.append ("-l")
	    if self.onesession:
		cmd.append ("-d")
		cmd.append (self.fd)

	    for a in args:
		if isinstance(a,list):
		    for x in a:
			cmd.append(x)
		else:
		    cmd.append(a)
	    if self.verbose:
		print "getINDI Sending", cmd
	    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
	    (props_string, errs) = p.communicate()
	    if (p.returncode):
		raise PyINDIException("getINDI exit code " + str(p.returncode) + ": " + errs)
	    if (errs):
		raise PyINDIException("getINDI errors: " + errs)

	    if showlabels:
		for l in sorted(props_string.splitlines()):
		    print l

	    else:
		props_dict = _crackGIP (props_string)
		if len(props_dict) == 1:
		    value = props_dict[props_dict.keys()[0]]
		    if self.verbose:
			print "getINDI Returning", value
		    return _giType(value)
		else:
		    for p in props_dict.keys():
			props_dict[p] = _giType(props_dict[p])
		    if self.verbose:
			print "getINDI Returning", props_dict
		    return props_dict



	# fetch an INDI BLOB, assummed to be a FITS file
	def getFITS (self, *args, **kwargs):
	    """Fetch an INDI BLOB, presumed to be a FITS file, as a pyfits object

	    Optional named arguments:
		timeout: max time to wait for FITS file, seconds, default 10

	    Usage:
		Call with one arg:
		    wait for and return the given BLOB property as a pyfits object.
		    Example: getFITS ("A.B.BLOB")
		Call with two args:
		    arg 1 is BLOB, arg 2 is trigger property:
			issue trigger then wait for BLOB
		    Example: getFITS ("A.B.BLOB", "A.X.Trigger")
	    """

	    # get desired timeout
	    to = str(kwargs.pop('timeout', 10))

	    if len(args) == 1:
		# just wait for blob synchronously
		blob = list(args)[0]
		cmd = ["getINDI", "-h", self.host, "-p", self.port, "-t", to, "-B"]
		cmd.append(blob)
		if self.onesession:
		    cmd.append ("-d")
		    cmd.append (self.fd)
		if self.verbose:
		    print "getFITS Sending", cmd
		p = Popen(cmd, stdout=PIPE, stderr=PIPE)
		(props_string, errs) = p.communicate()
		if (p.returncode):
		    raise PyINDIException("getINDI exit code " + str(p.returncode) + ": " + errs)
		if (errs):
		    raise PyINDIException("getINDI errors: " + errs)
		if self.verbose:
		    print "getFITS Returning", blob
		return (self._openFITS (blob))

	    elif len(args) == 2:
		# issue async get, issue set trigger, then collect BLOB
		al = list(args)
		blob = al[0]
		trigger = al[1]
		cmd = ["getINDI", "-h", self.host, "-p", self.port, "-t", to, "-B", blob]
		if self.onesession:
		    cmd.append ("-d")
		    cmd.append (self.fd)
		if self.verbose:
		    print "getFITS Sending", cmd
		p = Popen(cmd, stdout=PIPE, stderr=PIPE)
		self.setINDI (trigger, wait=True, timeout=to)
		(props_string, errs) = p.communicate()
		if (p.returncode):
                    raise PyINDIException("getINDI exit code " + str(p.returncode) + ": " + errs)
		if (errs):
		    raise PyINDIException("getFITS errors: " + errs)
		if self.verbose:
		    print "getFITS Returning", blob
		return (self._openFITS (blob))

	    else:
		raise PyINDIException ("getFITS with " + str(len(args)) + " args")


	# return dictionary of a set of related properties.
	def getKeys (self, device, prop, key, nunique=0):
	    """ Given an INDI property that is sent multiple times with different values for a given
	    element, collect each property into a dictionary indexed by each unique key value.

	    Required named arguments:
		device: INDI device name, no default
		prop:   property name to inspect, no default
		key:    element to use for dictionary key, no default

	    Options named arguments:
		nunique: skip any possible dups to insure exactly this many keys

	    Example:
		read each unique ID ROIStats from the NOMIC science camera:
		    rois = pi.getKeys (device='NOMIC', prop='ROIStats', key='ID', nunique=3)
		    for id in rois:
			print id, rois[id]['Sum']

	    """

	    # open pipe to getINDI
	    property = device + '.' + prop + '.*'
	    cmd = ['getINDI', '-h', self.host, '-p', self.port, '-m', '-t', '0', '-f', property]
	    if self.verbose:
		print "getKeys running ", cmd
	    try:
		g = Popen (cmd, stdout=PIPE)
	    except OSError:
		raise Exception ("Can not run " + str(cmd))

	    # read each line, building up rois{}, until find all
	    id = False
	    rois = {}
	    # for line in g.stdout:		# don't use this form, it does internal buffering
	    for line in iter(g.stdout.readline, b''):
		chomp = line.rstrip('\n')
		(p, v) = chomp.split('=')
		(d, n, e) = p.split('.')
		if self.verbose:
		    print "getKeys parsing ", chomp
		if e == key:
		    id = v
		    if nunique > 0 and len(rois) == nunique:
			break
		    if nunique == 0 and id in rois:
			break;
		    rois[id] = {}
		if id != False:
		    rois[id][e] = v
	    g.kill()
	    return rois


	# set INDI properties
	def setINDI (self, *args, **kwargs):
	    """Set one or more INDI properties, then optionally wait while Busy.

	    Optional named arguments:
		wait:    whether to wait or just return immediately, boolean, default True
		timeout: max time to wait if using wait, seconds, default 10

	    Usage:
		Call with any number or combinations of:
		    dictionary defining property name:value pair
		    separate name and value pair
		    atomic shorthand

		Examples:
		    setINDI("A.B.C", 1)                      # separate name, value pair
		    setINDI("A.B.C", 1, timeout=60)          # same, but with longer timeout
		    setINDI("A.B.C", 2, "X.Y.Z", 3)          # two separate pairs
		    setINDI({"A.B.C":2, "X.Y.Z":3})          # same but in a dictionary
		    setINDI({"M.N.O":2, "M.N.P":4})          # atomic setting in dictionary
		    setINDI("M.N.O=2;P=4")                   # same using shorthand
		    setINDI("X.Y.Filter", "Red", wait=False) # sync later with waitINDI()
	    """

	    # convert v to INDI string value
	    def _toIS (v):
		if isinstance(v,bool):
		    if v:
			v = "On"
		    else:
			v = "Off"
		elif not isinstance(v,str):
		    v = str(v)
		return v

	    # check whether to wait and how long
	    wait = kwargs.pop('wait', True)
	    to = str(kwargs.pop('timeout', 10))

	    # collect name,value pairs and/or dictionaries into nvp
	    # atomic settings go directly into specs
	    n = None
	    nvp = {}
	    specs = []
	    for a in args:
		if isinstance(a,dict):
		    nvp.update(a)
		else:
		    a = _toIS(a)
		    if a.find("=") >= 0:
			specs.append (a)
		    elif n == None:
			n = a
		    else:
			nvp[n] = a
			n = None
	    if n != None:
		raise PyINDIException("No value paired with " + n)

	    # add d to specs, combining any atomic groups
	    spec = ""
	    for n in sorted(nvp.keys()):
		n = _toIS(n)
		v = _toIS(nvp[n])
		comp = n.split(".")
		if spec == "":
		    spec = n + "=" + v
		elif prev[0] == comp[0] and prev[1] == comp[1]:
		    spec = spec + ";" + comp[2] + "=" + v
		else:
		    specs.append(spec)
		    spec = n + "=" + v
		prev = comp[:]
	    if spec != "":
		specs.append(spec)

	    # start the command line
	    cmd = ["setINDI", "-h", self.host, "-p", self.port, "-t", to]
	    if self.onesession:
		cmd.append ("-d")
		cmd.append (self.fd)

	    # add -w if waiting
	    if wait:
		cmd.append ("-w")

	    # add all specs to cmd
	    cmd += specs

	    # do it
	    if self.verbose:
		print "setINDI Sending", cmd
	    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
	    (props_string, errs) = p.communicate()
	    if (p.returncode):
		raise PyINDIException("setINDI exit code " + str(p.returncode) + ": " + errs)
	    if (errs):
		raise PyINDIException("setINDI errors: " + errs)

	    if self.verbose:
	        print "setINDI Returning"


	# wait while property is Busy
	def waitINDI (self, *args, **kwargs):
	    """Wait while INDI properties are Busy.

	    Optional named arguments:
		timeout: max time to wait, seconds, default 10

	    Usage:
		Call with one or more property name strings or lists of same:
		    wait until all property states are other than Busy

		Examples:
		    waitINDI ("A.B.C", "D.E.F")
		    waitINDI (["A.B.C", "D.E.F"])
	    """

	    def _doWait (a):

		# start the command line
		cmd = ["evalINDI", "-h", self.host, "-p", self.port, "-t", str(kwargs.pop('timeout', 10)), "-o", "-w"]
		if self.onesession:
		    cmd.append ("-d")
		    cmd.append (self.fd)

		# crack open the property, set element expression to _STATE != Busy
		p = a.split(".")
		cmd.append ('"' + p[0] + '.' + p[1] + '._STATE"!=2')
		if self.verbose:
		    print "waitINDI Sending", cmd

		# wait until STATE is not Busy
		p = Popen(cmd, stdout=PIPE, stderr=PIPE)
		(props_string, errs) = p.communicate()
                if (p.returncode):
                    raise PyINDIException("evalINDI exit code " + str(p.returncode) + ": " + errs)
		if (errs):
		    raise PyINDIException("waitINDI errors: " + errs)

		# check last result for anything other than 1
		l = props_string.splitlines()
		if l[len(l)-1].split("=")[1] != "1":
		    raise PyINDIException(a + " reports error")

	    for a in args:
		if isinstance(a, list):
		    for x in a:
			_doWait(x)
		else:
		    _doWait(a)

	    if self.verbose:
		print "waitINDI Returning"

	# wait until expression evaluates to true
	def evalINDI (self, *args, **kwargs):
	    """Wait until an expression of INDI properties evaluates as true.

	    Optional named arguments:
		timeout: max time to wait, seconds, default 10

	    Usage:
		Call with a string containing any arithmetic expression using INDI properties
		in double-quotes:
		    wait until expression evaluates as true or times out.

		Example:
		    evalINDI ('"A.B.C" < 5');
	    """

	    def _doEval (a):

		# start the command line
		cmd = ["evalINDI", "-h", self.host, "-p", self.port, "-w", "-t", str(kwargs.pop('timeout', 10))]
		if self.onesession:
		    cmd.append ("-d")
		    cmd.append (self.fd)

		# add the expression
		cmd.append (a)
		if self.verbose:
		    print "evalINDI Sending", cmd

		# wait for evalINDI to exit when the expression is true or timed out
		p = Popen(cmd, stdout=PIPE, stderr=PIPE)
		(props_string, errs) = p.communicate()
                if (p.returncode):
                    raise PyINDIException("evalINDI exit code " + str(p.returncode) + ": " + errs)
		if (errs):
		    raise PyINDIException("evalINDI errors: " + errs)

	    for a in args:
		if isinstance(a, list):
		    for x in a:
			_doEval(x)
		else:
		    _doEval(a)

	    if self.verbose:
		print "evalINDI Returning"

	# handy dictionary pretty-print
	def ppD (self, d, pattern="*"):
	    """Handy Pretty-Printer for any Dictionary

	    Optional second argument:
	        pattern: file name glob pattern for dictionary key, default "*"

	    Usage:
		Call with any dictionary, see each entry sorted in the form n=v, one per line

		Examples:
		    ppD(dict)			# print all entries in dict
		    ppD(dict, "Time.*.*")	# print only entries whose key matchs the pattern
	    """

	    # scan once to get longest key
	    maxl = 0
	    for k in d:
		if fnmatch.fnmatch(k,pattern) and len(k) > maxl:
		    maxl = len(k)

	    # scan again to print in sorted order
	    for k in sorted(d):
		if (fnmatch.fnmatch(k,pattern)):
		    print "%-*s = %s" % (maxl, k, d[k])

	# handy image centroider
	def centroid (self, pixels):
	    """Handy method to compute centroid of an image.

	    Usage:
		Call with the data from a FITS image, returns tuple with
		    (centroid_row, centroid_col)

		Example:
		    from pyindi import *
		    import numpy
		    pi = PyINDI("lbti-lmircam")
		    f = pi.getFITS("LMIRCam.DisplayImage.File", "LMIRCam.GetDisplayImage.Now", "On")
		    pixels = f[0].data
		    (r, c) = pi.centroid(pixels)
		    v = pixels[r,c]
		    sigma = (v - pixels.mean())/pixels.std()
	    """

	    # get shape and sum
	    (rows, cols) = pixels.shape
	    sum = pixels.astype(numpy.float64).sum()	# float32 can easily underflow

	    # find centroids
	    r = pixels.transpose().dot(numpy.array(range(rows))).sum()/sum
	    c = pixels.dot(numpy.array(range(cols))).sum()/sum
	    return (int(r+0.5), int(c+0.5))

	# handy region cut
	def cutRegion (self, pixels, r, c, w, h):
	    """Handy method to extract a region from an array of pixels

	    Usage:
		Call with a 2D ndarray, starting row, col and width, height:
		    return the given region as a new image

	    Example:

		# compute centroid of a narrow 50x200 region starting at 10,20:
		roi = pi.cutRegion (pixels, 10, 20, 50, 200)
		(r, c) = pi.centroid(roi)
		v = roi[r,c]
	    """

	    return (pixels[r:r+h, c:c+w])

	# handy H:M:S
	def hms (self, h):
	    if (h < 0):
		isneg = True
		h = -h
	    else:
		isneg = False
	    hi = int(h)
	    m = 60*(h-hi)
	    mi = int(m)
	    s = 60*(m-mi)
	    r = "%d:%02d:%06.3f" % (hi, mi, s)
	    if (isneg):
		r = "-" + r
	    return r


# sample unit test if run alone
if __name__ == '__main__':
    # pi = PyINDI()
    # pi = PyINDI(host='lbti-web', onesession=True)
    pi = PyINDI(host='lbti-nomic')
    # pi = PyINDI("127.0.0.1")

    # help (pi)

    for i in range(10):
	rois = pi.getKeys (device='NOMIC', prop='ROIStats', key='ID', nunique=3)
	print
	for id in rois:
	    print id, rois[id]['Sum']

    # get exactly one
    # print pi.getINDI('Time.Now.UTC')

    # get all, print a subset
    # pd = pi.getINDI(timeout=1)
    # print pd['Time.Now.UTC'], pd['Time.Now.LST']

    # get a subset with wild card
    # pd = pi.getINDI('Time.Now.*')
    # print pd['Time.Now.UTC'], pd['Time.Now.LST']

