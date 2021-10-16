#! /usr/bin/env python3
# GPFPlot Library for Parsing Input File
# Last modified: 2020-07-09

MORB = ["ORB","HFORB"]
DENS = ["QC","INT","TOT"]
MISC = ["PROMPT"]
MODES= MORB + DENS + MISC

def getLine(source, string, pos=0, none=True, comments=True):
	"""Get the line index of the ($pos+1)-th
       occurrence of $string in $source."""
	count=0
	for index in range(len(source)):
		if string in source[index]:
			if source[index][0] == "#":
				continue # Ignore comments
			if pos == count:
				return index
			count += 1
	if none:
		return None
	else:
		raise ValueError('string not found')

def truefalse(s):
  return s.lower() in ("yes", "true", "t", "y", "1")

def parse(file, string, default=""):
	try:
		# string+"=" prevents finding variables instead of values 
		line = getLine(file, string+"=")
		result = file[line].split("=")[1]
	except:
		result = default
	return result	

def parse1(infile):
	out_file  = parse(infile,  "out_file")
	dens_file = parse(infile, "dens_file")
	mode = parse(infile, "mode", "PROMPT")
	return out_file, dens_file, mode

def parse2(infile):
	plane = parse(infile,  "plane", "xy")
	gridp = int( parse(infile, "gridp", "50"))
	Xlim = tuple( map(float, parse(infile, "xlim", "-1,1").split(",")))
	Ylim = tuple( map(float, parse(infile, "ylim", "-1,1").split(",")))
	offset = float( parse(infile, "offset", "0"))
	return plane, gridp, Xlim, Ylim, offset

def parse3(infile, mode):
	c_fill  = truefalse( parse(infile, "c_fill", "t"))
	if mode in DENS:
		color_f = "cm." + parse(infile, "color", "seismic_r")
		min_c = float( parse(infile, "min_c", "-.1")) 
		max_c = float( parse(infile, "max_c",  ".1")) 
	else:
		color_f = "cm." + parse(infile, "color_f", "coolwarm")
		min_c = float( parse(infile, "min_c", "-1")) 
		max_c = float( parse(infile, "max_c",  "1")) 
	color_l = parse(infile, "color_l")
	c_lines = truefalse( parse(infile, "c_lines", "t"))
	c_label = truefalse( parse(infile, "c_label", "false"))
	auto_c= truefalse( parse(infile, "auto_c", "false"))
	nconts= int( parse(infile, "nconts",  "10"))
	p_unit= parse(infile, "p_unit", "angs")
	draw_atom = parse(infile, "draw_atom","plane")
	draw_name = truefalse( parse(infile, "draw_name", "t"))
	return c_fill, color_f, c_lines, color_l, c_label,\
           min_c, max_c, nconts, auto_c,\
           p_unit, draw_atom, draw_name

def parse4(infile):
	save_list = parse(infile,  "save_png").split(" ")
	if len(save_list) == 0:
		save_png = ""
	elif len(save_list) >= 1:
		save_png = save_list[0]
	if len(save_list) > 1:
		dpipng = int( save_list[1] )
	else:
		dpipng = 300
	save_list = parse(infile,  "save_eps").split(" ")
	if len(save_list) == 0:
		save_eps = ""
	elif len(save_list) >= 1:
		save_eps = save_list[0]
	if len(save_list) > 1:
		dpieps = int( save_list[1] )
	else:
		dpieps = 300
	save_txt = parse(infile, "save_txt")
	return save_png, save_eps, dpipng, dpieps, save_txt
