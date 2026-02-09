import sys
from os import makedirs
from pathlib import Path
import numpy as np
import pandas as pd
from bokeh.models import ColumnDataSource, Whisker, Range1d, Label
from bokeh.plotting import figure, save
from bokeh.layouts import gridplot
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline, PchipInterpolator, Akima1DInterpolator
from scipy.interpolate import make_splrep
from sklearn import linear_model
# installed as "scikit-learn"
# from bokeh.io import export_svg
from bokeh.io import export_png
import regex as re
re.DEFAULT_VERSION = re.VERSION1
# Needs pip install regex. normal re doesn't support varying length of look-behind patterns.
import datetime

drugcolor = ("#4477AA", "#228833", "#CCBB44", "#EE6677")
# The colors "#4477AA", "#228833", "#CCBB44", "#EE6677" are from https://personal.sron.nl/~pault/
drugdash = ("solid", "solid", "solid", "solid")
# Line pattern, such as "solid", "dashed", "dotted", or "dotdash".
lwidth = 4
# line width
wlwidth = 2.5
# whisker line width
walpha = 0.6
# whisker alpha
addtitle = True
# True or False, not in quotes
# If true, file name minus filetags will be used as title.
showtoolbar = False
# Show or hide the Bokeh web toolbar.
reg_type = "pchip"
# How to model the viability data: "none", "pchip", "akima", "4PL", "5PL", "b-spline", or "catmull_rom" (with the quotes)
# If "none" is selected, will use connected line segments to calculate IC50s.
# Only "pchip" is known to be functional as of now.

def mttg(datafile):
    # datafile = sys.argv[1]
    # need to catch error when the file is already open in excel
    # datafile = filetemp
    datafilename = Path(datafile).stem
    print(datafilename)

    filetags = re.findall(r"(?<=(?:(?:\.\.mb\+)|(?:(?:\.\.mb\+)(?:[^(\+)(\.\.)]*\+)+)))[^(?:\+)(?:\.\.)]*", datafilename)
    filetags
    # This gives a list of the tags in the file name. For example,
    # in a filename that has this string somewhere in it, 
    #   "..mb+NC=H04H06+PC=H10H12+PB=H10H12+PE=H10H12..",
    # it will return ['NC=H04H06', 'PC=H10H12', 'PB=H10H12', 'PE=H10H12'].
    # Andvanced Renamer and similar tools are useful for adding tags.

    if not list(filter(lambda v: re.match(r"^TF=", v), filetags)):
        # lambda and filter() usage as in https://stackoverflow.com/a/64143128.
        # The "not" checks if list is empty (in case of using default template) as in https://stackoverflow.com/questions/53513/how-do-i-check-if-a-list-is-empty
        templatefile = str(sys.path[0]) + "/template.xlsx"
    else:
        templatefile = str(sys.path[0]) + "/template" + str(list(filter(lambda v: re.match(r"^TF=", v), filetags))[0][3]) + ".xlsx"

    # Currently, template must have a single first column of concentrations
    # and the full 8x12 plate represented to its right. Empty wells may
    # be left empty currently.

    tp = pd.read_excel(
        io=templatefile,
        usecols="B:N",
        skiprows=37,
        # need to make this automatically adjustable based on loc of "conc" and maybe "end"
        header=None
    )

    md = pd.read_excel(
        io=datafile,
        usecols="C:N",
        skiprows=37,
        header=None
    )


    tfac = {"row": int(tp.shape[0] / 8), "col": int((tp.shape[1] - 1) / 12)}
    # 1 is subtracted in "col" for concs column
    def tfplate(plate_positions):
        return range((ord(plate_positions[0]) - 65) * tfac["row"], ((ord(plate_positions[3]) - 65) * tfac["row"]) + (tfac["row"] - 1) + 1), range((int(plate_positions[1:3]) - 1) * tfac["col"] + 1, ((int(plate_positions[4:6]) - 1) * tfac["col"]) + (tfac["col"] - 1 + 1 + 1))
    # Takes a plate position of the form "X##X##" (including, e.g., A10A10 for one well) for use with tp.iloc[var with tfplate output].

    def mplate(plate_positions):
        return range((ord(plate_positions[0]) - 65) * tfac["row"], ((ord(plate_positions[3]) - 65) * tfac["row"]) + (tfac["row"] - 1) + 1), range((int(plate_positions[1:3]) - 1) * tfac["col"], ((int(plate_positions[4:6]) - 1) * tfac["col"]) + (tfac["col"] - 1 + 1))
    # Takes a plate position of the form "X##X##" (including, e.g., A10A10 for one well) for use with tp.iloc[var with tfplate output].


    tMCold = np.where(tp == "VC.pos"), np.where(tp == "VC.neg")
    print(tMCold)
    # np.where usage like in https://stackoverflow.com/a/49669857
    # tMCold = range(min([min(tMCold[0][0]), min(tMCold[1][0])]), max([max(tMCold[0][0]), max(tMCold[1][0])]) + 1), range(min([min(tMCold[0][1]), min(tMCold[1][1])]), max([max(tMCold[0][1]), max(tMCold[1][1])]) + 1)
    # Guestimate the position of the template file's control wells. With hope it has the same dimensions and orientation as that specified by the file tag.
    if list(filter(lambda v: re.match(r"^MC=", v), filetags)):
        tMCtrl = tfplate(list(filter(lambda v: re.match(r"^MC=", v), filetags))[0][3:])
        # lambda and filter() usage as in https://stackoverflow.com/a/64143128.
        # M for mortality.
        # tp = tp.replace(["MC.pos", "MC.neg"], np.nan)
        # tp.iloc[tMCtrl] = tp.iloc[tMCold]
        # tp.iloc[tMCold] = np.nan
        # Then, move those values, as in https://stackoverflow.com/a/70020841
        mtMCtrl = mplate(list(filter(lambda v: re.match(r"^MC=", v), filetags))[0][3:])
        mort_avg = np.nanmean(md.iloc[mtMCtrl].iloc[0,:] - md.iloc[mtMCtrl].iloc[1,:])
        # temporaryyy
    else:
        mort_avg = np.nanmean(md.iloc[tMCold].iloc[0,:] - md.iloc[tMCold].iloc[1,:])
    
    tVCold = np.where(tp == "VC.pos"), np.where(tp == "VC.neg")
    # np.where usage like in https://stackoverflow.com/a/49669857
    # tVCold = range(min([min(tVCold[0][0]), min(tVCold[1][0])]), max([max(tVCold[0][0]), max(tVCold[1][0])]) + 1), range(min([min(tVCold[0][1]), min(tVCold[1][1])]), max([max(tVCold[0][1]), max(tVCold[1][1])]) + 1)
    # Guestimate the position of the template file's control wells. With hope it has the same dimensions and orientation as that specified by the file tag.
    if list(filter(lambda v: re.match(r"^VC=", v), filetags)):
        # V for viability.
        tVCtrl = tfplate(list(filter(lambda v: re.match(r"^VC=", v), filetags))[0][3:])
        # tp.iloc[tVCtrl] = tp.iloc[tVCold]
        # tp.iloc[tVCold] = np.nan
        # Then, move those values, as in https://stackoverflow.com/a/70020841
        mtVCtrl = mplate(list(filter(lambda v: re.match(r"^VC=", v), filetags))[0][3:])
        viab_avg = np.nanmean(md.iloc[mtVCtrl].iloc[0,:] - md.iloc[mtVCtrl].iloc[1,:])
        # temporaryyy
    else:
        viab_avg = np.nanmean(md.iloc[tVCold].iloc[0,:] - md.iloc[tVCold].iloc[1,:])

    if list(filter(lambda v: re.match(r"^RM\d=", v), filetags)):
        # if the list of rmv tags (wells in X##X## form to be removed) is not empty, return true
        rmvs_wellform = list(filter(lambda v: re.match(r"^RM\d=", v), filetags))
        rmvs = []
        for i in range(len(rmvs_wellform)):
            rmvs.append(tfplate(rmvs_wellform[i][4:]))
            tp.iloc[rmvs[i]] = np.nan

    # Get the drug names with Regex
    drugphrase = re.search(r"(?<= )[^ ]+(?= \d\dhrs)", datafilename).group(0)
    drugname = re.findall(r"\d{4}|\w{3}", str(drugphrase))
    numS = len(drugname)


    print(tp)
    concslist = tp.iloc[:,0].dropna().drop_duplicates()
    tp_noconc = tp.drop(tp.columns[0], axis=1).reset_index(drop = True)
    print(tp_noconc)
    # list of data frames, one for each drug's sets of abs readings
    templsampls = tp_noconc.stack().unique()
    # get list of all unique template cell names
    sdata = dict()
    numreps = dict()
    tsPosList = list(filter(lambda v: re.match(r"\d\.pos", v), filter(lambda b: b==b, templsampls)))
    # nans are filtered out, per https://stackoverflow.com/a/29679784
    print(tsPosList)
    for i in tsPosList:
        # data must have "pos" ABS values, but "neg" ones are optional
        # But this does not yet handle multiple "pos" readings for one well
        numreps[i] = max(np.where(tp == i, 1, 0).sum(axis=1))
        # "where" usage as in https://stackoverflow.com/a/73557704
        # this gets the maximum number of replicates to provide columns for
        sdata[i] = pd.DataFrame(np.nan, index = list(range(len(concslist))), columns = list(range(numreps[i] + 1)))
        sdata[i].rename(columns = {sdata[i].columns[numreps[i]]: "avg"}, inplace = True)
        # empty creation as in https://stackoverflow.com/a/30053507
        # sdata[i] = sdata[i].reindex(columns = list(range(numreps[i])))
        # sdata[i].index = range(len(concslist))
        for j in range(len(concslist)):
            # sdata[i].loc[j, "mk"] = dict()
            tempincr = 0
            for u in range(len(tp_noconc.iloc[j * 2, :])):
                if tp_noconc.iloc[j * 2, u] == i:
                    sdata[i].iloc[j, tempincr] = ((md.iloc[j * 2, u] - md.iloc[(j * 2) + 1, u]) - mort_avg) / (viab_avg - mort_avg)
                    # the neglection of .neg managment is temporarryyyy
                    tempincr += 1
            sdata[i].loc[j, "avg"] = np.nanmean(sdata[i].iloc[j, 0:(numreps[i])])

    for i in range(numS):
        sdata[i] = sdata.pop(str(i + 1) + ".pos")
        numreps[i] = numreps.pop(str(i + 1) + ".pos")

    # Whiskers
    def addwhbounds(samplenum):
        vstd = sdata[samplenum].iloc[:, :-1].std(axis = 1, ddof = 1)
        # get all columns except the last one (and get std dev)
        lower = sdata[samplenum]["avg"] - vstd
        higher = sdata[samplenum]["avg"] + vstd
        return [lower, higher]

    for i in range(numS):
        sdata[i]["lowerW"], sdata[i]["upperW"] = addwhbounds(i)

    def spline_4p( t, p_1, p0, p1, p2 ):
        """ Catmull-Rom
            (Ps can be numpy vectors or arrays too: colors, curves ...)
        """
            # wikipedia Catmull-Rom -> Cubic_Hermite_spline
            # 0 -> p0,  1 -> p1,  1/2 -> (- p_1 + 9 p0 + 9 p1 - p2) / 16
        # assert 0 <= t <= 1
        return (
            t*((2-t)*t - 1)   * p_1
            + (t*t*(3*t - 5) + 2) * p0
            + t*((4 - 3*t)*t + 1) * p1
            + (t-1)*t*t         * p2 ) / 2
        # this function from https://stackoverflow.com/a/1295081
    def log4pl(x, A, B, C, D):
        return(((A-D)/(1.0+((x/C)**B))) + D)
    def log5pl(x, A, B, C, D, E):
        return((A-D)/((1.0+((x/C)**B))**E) + D)
        # above 4 lines from https://medium.com/@tentotheminus9/elisa-analysis-in-python-deb8c6ed91db
    # def residuals(p, y, x):
    #     """Deviations of data from fitted 4PL curve"""
    #     A,B,C,D = p
    #     err = y-log4pl(x, A, B, C, D)
    #     return err
    if reg_type == "b-spline" or reg_type == "pchip" or reg_type == "akima":
        # need to try https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-for-a-dataset instead
        # regxvals = pd.concat([mdv1["conc"] - 0.00000001, mdv1["conc"], mdv1["conc"] + 0.00000001])
        regxvals = concslist.reset_index(drop = True)
        print(regxvals)
        # regxvals = pd.concat([mdv1["conc"] - 3, mdv1["conc"], mdv1["conc"] + 3])
        # regxvals = pd.concat([mdv1["conc"], mdv1["conc"], mdv1["conc"]])
    else:
        # regxvals = pd.concat([mdv1["conc"], mdv1["conc"], mdv1["conc"]])
        regxvals = np.repeat(concslist, 3)
    regvals = [0] * numS
    regparams = pd.DataFrame()
    # regressions = pd.DataFrame(np.logspace(-0.85387196432, 2, 1000), columns = ["concs"])
    regressions = [0] * numS
    recovbyreg = pd.DataFrame()
    # avgrelativefit = [None] * numS
    avgregressions = [0] * numS
    inco50est = [0] * numS
    inco50estest = [0] * numS
    inco50 = [0] * numS
    if reg_type == "b-spline":
        b_splines = [0] * numS
    for i in range(numS):
        regvals[i] = [0] * numreps[i]
        regressions[i] = [0] * (numreps[i] + 2)
        # regvals = [0] * numreps[i]
        for j in range(numreps[i]):
            if reg_type == "b-spline" or reg_type == "pchip" or reg_type == "akima":
                regvals[i][j] = pd.concat([sdata[i][j].rename("vals"), regxvals.rename("concs")], axis = 1)
            else:
                regvals[i] = pd.concat([pd.concat([mdv1.iloc[:, i * 4], mdv1.iloc[:, i * 4 + 1], mdv1.iloc[:, i * 4 + 2]]).rename("vals"), regxvals.rename("concs")], axis = 1)
            regvals[i][j] = regvals[i][j][regvals[i][j]["vals"].notna()]
            if reg_type == "4PL":
                regparams[("Drug " + str(i) + " Params")], _ = curve_fit(log4pl, regvals[i]["concs"].to_numpy(), regvals[i]["vals"].to_numpy(), maxfev = 50000)
                    # above line from https://medium.com/@tentotheminus9/elisa-analysis-in-python-deb8c6ed91db
                # regparams[("Drug " + str(i) + " Params")] = leastsq(residuals, p0, args = (regvals[i]["vals"].to_numpy(), regvals[i]["concs"].to_numpy()))[0]
                regressions[i] = pd.DataFrame(np.logspace(np.log10(regvals[i]["concs"].min()), 2, 1000), columns = ["concs"])
                regressions[i]["regvals"] = log4pl(regressions[i]["concs"], regparams.iloc[0,i], regparams.iloc[1,i], regparams.iloc[2,i], regparams.iloc[3,i])
                regvals[i]["recov" + str(i)] = log4pl(regvals[i]["concs"], regparams.iloc[0,i], regparams.iloc[1,i], regparams.iloc[2,i], regparams.iloc[3,i])
                # print((regressions[i]["regvals"] - 0.5).abs().min)
                # print(regressions[i].loc[(regressions[i]["regvals"] - 50).abs().idxmin()])
                # print((regressions[i]["regvals"] - 0.5).abs().idxmin())
                # Part of this line from https://stackoverflow.com/a/60116914
                # inco50[i] = regparams.iloc[2,i]*((50 - regparams.iloc[0,i])/(regparams.iloc[3,i] - 0.5))**(1/regparams.iloc[1,i])
                # inverse according to Wolfram Alpha
            elif reg_type == "5PL":
                regparams[("Drug " + str(i) + " Params")], _ = curve_fit(log5pl, regvals[i]["concs"].to_numpy(), regvals[i]["vals"].to_numpy(), maxfev = 50000)
                    # above line from https://medium.com/@tentotheminus9/elisa-analysis-in-python-deb8c6ed91db
                # regparams[("Drug " + str(i) + " Params")] = leastsq(residuals, p0, args = (regvals[i]["vals"].to_numpy(), regvals[i]["concs"].to_numpy()))[0]
                regressions[i] = pd.DataFrame(np.logspace(np.log10(regvals[i]["concs"].min()), 2, 1000), columns = ["concs"])
                regressions[i]["regvals"] = log5pl(regressions[i]["concs"], regparams.iloc[0,i], regparams.iloc[1,i], regparams.iloc[2,i], regparams.iloc[3,i], regparams.iloc[4,i])
                regvals[i]["recov" + str(i)] = log5pl(regvals[i]["concs"], regparams.iloc[0,i], regparams.iloc[1,i], regparams.iloc[2,i], regparams.iloc[3,i], regparams.iloc[4,i])
                # Part of this line from https://stackoverflow.com/a/60116914
                # inco50[i] = regparams.iloc[2,i] * (((regparams.iloc[3,i] - regparams.iloc[0,i])/(regparams.iloc[3,i] - 50))**(1/regparams.iloc[4,i]) - 1)**(1/regparams.iloc[1,i])
                # inverse according to Wolfram Alpha
            elif reg_type == "catmull_rom":
                regparams[("Drug " + str(i) + " Params")], _ = curve_fit(spline_4p, regvals[i]["concs"].to_numpy(), regvals[i]["vals"].to_numpy(), maxfev = 50000)
                    # above line from https://medium.com/@tentotheminus9/elisa-analysis-in-python-deb8c6ed91db
                regressions[i] = pd.DataFrame(np.logspace(np.log10(regvals[i]["concs"].min()), 2, 1000), columns = ["concs"])
                regressions[i]["regvals"] = spline_4p(regressions[i]["concs"], regparams.iloc[0,i], regparams.iloc[1,i], regparams.iloc[2,i], regparams.iloc[3,i])
                regvals[i]["recov" + str(i)] = spline_4p(regvals[i]["concs"], regparams.iloc[0,i], regparams.iloc[1,i], regparams.iloc[2,i], regparams.iloc[3,i])
            elif reg_type == "pchip":
                regressions[i][j] = pd.DataFrame(np.logspace(np.log10(regvals[i][j]["concs"].min()), 2, 1000), columns = ["concs"])
                print(regvals[i][j])
                monospline = PchipInterpolator(regvals[i][j].sort_values("concs")["concs"].to_numpy(), regvals[i][j].sort_values("concs")["vals"].to_numpy())
                regressions[i][j]["regvals"] = monospline(regressions[i][j]["concs"])
                regvals[i][j]["recov" + str(i) + str(j)] = monospline(regvals[i][j]["concs"])
            elif reg_type == "akima":
                regressions[i] = pd.DataFrame(np.logspace(np.log10(regvals[i]["concs"].min()), 2, 1000), columns = ["concs"])
                monospline = Akima1DInterpolator(regvals[i].sort_values("concs")["concs"].to_numpy(), regvals[i].sort_values("concs")["vals"].to_numpy())
                regressions[i]["regvals"] = monospline(regressions[i]["concs"])
                regvals[i]["recov" + str(i)] = monospline(regvals[i]["concs"])
            elif reg_type == "b-spline":
                # b_splines[i] = make_smoothing_spline(regvals[i].sort_values("concs")["concs"].to_numpy(), regvals[i].sort_values("concs")["vals"].to_numpy())
                # The spline is cubic as it has k=3 (degree of 3) by default.
                b_splines[i] = make_splrep(regvals[i].sort_values("concs")["concs"].to_numpy(), regvals[i].sort_values("concs")["vals"].to_numpy(), k = 2, s = 100 * len(regvals[i].sort_values("concs")["concs"].to_numpy()))
                regressions[i] = pd.DataFrame(np.logspace(np.log10(regvals[i]["concs"].min()), 2, 1000), columns = ["concs"])
                # regressions[i]["regvals"] = sm.nonparametric.lowess(regvals[i]["vals"].to_numpy(), regvals[i]["concs"].to_numpy(), xvals = np.logspace(np.log10(regvals[i]["concs"].min()), 2, 1000))
                regressions[i]["regvals"] = b_splines[i](regressions[i]["concs"])
                regvals[i]["recov" + str(i)] = b_splines[i](regvals[i]["concs"])
                # regvals[i]["recov" + str(i)] = sm.nonparametric.lowess(regvals[i]["vals"].to_numpy(), regvals[i]["concs"].to_numpy(), return_sorted = False)
                # Part of this line from https://stackoverflow.com/a/60116914
            elif reg_type == "none":
                clf = linear_model.LinearRegression()
                # https://docs.scipy.org/doc/scipy/tutorial/interpolate/1D.html
                # noot finished yetttt
        regressions[i][-2] = pd.concat([regressions[i][x]["concs"] for x in range(numreps[i])], axis = 1)
        regressions[i][-1] = pd.concat([regressions[i][x]["regvals"] for x in range(numreps[i])], axis = 1)
        print(regressions[i][-1])
        avgregressions[i] = pd.DataFrame({"concs": regressions[i][-2].mean(axis = 1), "regvals": regressions[i][-1].mean(axis = 1)})
        # avgrelativefit[i] = np.average(regvals[i]["recov" + str(i)] / regvals[i]["vals"])
        # print(regressions[i].loc[regressions[i]["concs"] == 100, "regvals"].values[0])
        if avgregressions[i].loc[avgregressions[i]["concs"] == 100, "regvals"].values[0] <= 0.5:
            inco50estest[i] = [0] * numreps[i]
            for j in range(numreps[i]):
                inco50estest[i][j] = regressions[i][j].loc[(regressions[i][j]["regvals"] - 0.5).abs().idxmin()]["concs"]
            inco50est[i] = np.format_float_positional(avgregressions[i].loc[(avgregressions[i]["regvals"] - 0.5).abs().idxmin()]["concs"], precision = 4, unique = False, fractional = False, trim = "k")
            # lower number of digits produced (https://stackoverflow.com/a/58491045)
            inco50[i] = str(inco50est[i]) + " ± " + str(np.format_float_positional(np.nanstd(inco50estest[i], ddof = 1), precision = 3, unique = False, fractional = False, trim = "k"))
        else:
            inco50[i] = np.nan

    # print("avgfit with " + reg_type + " :")
    # print(avgrelativefit)
    print(inco50)
    # mcds = ColumnDataSource(data=mdv1)

    p = figure(x_axis_label="Drug concentration, μM", y_axis_label="Viability", x_axis_type="log")

    # Axis properties
    p.xaxis.ticker = concslist
    # Custom x axis tick locations
    # p.x_range = Range1d(xrangemin, xrangemax)
    p.axis.major_tick_in = 0
    p.axis.minor_tick_in = 0


    # print(mdv1[["t1avg", "t2avg", "t3avg", "t4avg"]].min().min())
    if min([min(sdata[x]["avg"]) for x in range(numS)]) <= 0:
        p.line([0.14, 100], [0, 0], line_width=lwidth / 2, color="black")
    # line at zero if lowest point is below zero
    p.line([0.14, 100], [0.5, 0.5], line_width=lwidth / 2, color="gray")
    # line at 50% viability
    via50label = Label(text_font_size="14px", text_color="gray", x=0.145, y=0.5025, x_units="data", y_units="data", text="50% viability", border_line_alpha=0, background_fill_alpha=0)
    p.add_layout(via50label)

    for i in range(numS):
        # p.line(x="conc", y=("t"+str(i+1)+"avg"), legend_label=(drugname[i] + " [IC50: " + str(inco50[i])[:4] + "]"), color=drugcolor[i], line_dash=drugdash[i], line_width=lwidth, level = "glyph", source = mcds)
        p.line(x=avgregressions[i]["concs"], y=avgregressions[i]["regvals"], legend_label=(drugname[i] + " [IC50: " + str(inco50[i]) + "]"), color=drugcolor[i], line_width=lwidth, level = "glyph", line_dash=drugdash[i])
        p.scatter(x=concslist, y=sdata[i]["avg"], line_color=drugcolor[i], line_width=2.5, color="white", level="overlay", size=7)
        # p.scatter(x="conc", y=("t"+str(i+1)+"avg"), line_color=drugcolor[i], line_width=2.5, color="white", level="overlay", size=7, source = mcds)


    # Legend properties
    p.legend.location = "bottom_left"
    p.legend.title = "Drugs"


    error = [None] * numS
    ersrc = [None] * numS
    print(sdata[0])
    for i in range(numS):
        # tmpx = [100.0, 33.0, 11.0, 3.7, 1.3, 0.4]
        ersrc[i] = ColumnDataSource(data=dict(base=concslist, upper=sdata[i]["upperW"], lower=sdata[i]["lowerW"]))
        # for j in range(len(sdata[i]["upperW"])):
        #     if isinstance(sdata[i].loc[j, "upperW"], float):
        #         tmpx.append(list(concslist)[j])
        error[i] = Whisker(base = "base", upper = "upper", lower = "lower", level = "annotation", line_width = wlwidth, source=ersrc[i])
        # at the same render levels, they are drawn in the order that they are added. otherwise: https://docs.bokeh.org/en/2.4.1/docs/user_guide/styling.html#setting-render-levels
        error[i].upper_head.line_alpha = walpha
        error[i].lower_head.line_alpha = walpha
        error[i].upper_head.line_width = wlwidth
        error[i].lower_head.line_width = wlwidth
        error[i].line_color = drugcolor[i]
        error[i].line_alpha = walpha
        error[i].upper_head.line_color = drugcolor[i]
        error[i].lower_head.line_color = drugcolor[i]
        error[i].upper_head.size = 17
        error[i].lower_head.size = 17
        p.add_layout(error[i])

    if addtitle == True:
        p.title=re.search(r".+?(?=\.\.mb)", datafilename).group(0)
        # p.title.text_font=textfont
        p.title.text_color='black'
        # p.title.text_font_size=titlefontsize_
        p.title.align='center'
    
    outmtt = [p, re.search(r".+?(?=\.\.mb)", datafilename).group(0), drugname, inco50]
    return outmtt


plots = []
comptable = pd.DataFrame()
for raw in Path(str(sys.path[0]) + r"\rawdata").glob('*.xlsx'):
    # from https://stackoverflow.com/a/43669828
    outmtt = mttg(raw)
    incotmpdict = {}
    incotmpdict["time (h)"] = re.search(r"\d\d(?=hrs)", outmtt[1]).group(0)
    for i in range(len(outmtt[2])):
        incotmpdict[outmtt[2][i]] = outmtt[3][i]
    incotmpdict
    incotmp = pd.DataFrame(incotmpdict, index=[0])
    comptable = pd.concat([comptable, incotmp], axis=0)
    plots.append(outmtt[0])

print(comptable)
rundate = str(datetime.datetime.now()).replace(":", ",")
makedirs("./output/" + rundate)
# print("/output/" + rundate + "/IC50s.csv")
comptable.to_csv("./output/" + rundate + "/IC50s.csv")

grid = gridplot(plots, ncols=3, width=675, height=675)
if showtoolbar == False:
    grid.toolbar_location = None
save(grid, filename=("./output/" + rundate + "/web.html"))
export_png(grid, scale_factor=2, filename=("./output/" + rundate + "/image.png"))

