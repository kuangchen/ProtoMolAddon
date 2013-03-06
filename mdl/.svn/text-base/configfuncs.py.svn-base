import sys

def check_envvar( name, die = False ):
    if os.environ.has_key(name):
        return os.environ[name]
    else:
        if die:
            print "Environment Variable Missing: " + name
            sys.exit(1)
        else:
            return None

def check_header( name, die = False ):
	if not conf.CheckCXXHeader( name ):
		if die:
			print "Header Missing: " + name
			sys.exit(1)
		else:
			return False
	else:
		return True

def check_library(conf,  name, die = False ):
	if not conf.CheckLib( name ):
		if die:
			print "Library Missing: " + name
			sys.exit(1)
		else:
			return False
	else:
		return True
