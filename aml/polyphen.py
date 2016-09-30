import os
import re
import pandas as pd
import urllib2
import cookielib
from matplotlib import pyplot
import matplotlib as mpl
from matplotlib.colors import ListedColormap

script_dir = os.path.dirname(os.path.abspath(__file__))


def create_polyphen_input(chrom, pos, ref, alt):
    """Returns a file with variant in format (e.g. chr1:1234567 A/G) required by PolyPhen-2 API.

    """
    input_file = "%s_%s.txt" % (chrom, pos)
    with open(input_file, "w+") as i:
        i.write("%s:%s %s/%s\n" % (chrom, pos, ref, alt))
    i.close()

    return input_file


def run_polyphen(input_file):
    """Connects to PolyPhen-2 using curl and outputs session information into a .txt file.

    Note: HumDiv and HumVar analysed separately.

    """
    os.system("curl "
              "-F _ggi_project=PPHWeb2 "
              "-F _ggi_origin=query "
              "-F _ggi_target_pipeline=1 "
              "-F MODELNAME=HumDiv "  # HUMDIV or HUMVAR
              "-F UCSCDB=hg19 "  # genome
              "-F SNPFUNC=m "  # functional SNP category, c or m
              "-F _ggi_batch_file=@%s "  # filename of variants, format "e.g. chr1:1158631 A/C/G/T"
              "-D - http://genetics.bwh.harvard.edu/cgi-bin/ggi/ggi2.cgi "
              "-c %s.humdiv.cookies.txt "
              "> %s.humdiv.txt" % (input_file, input_file[:-4], input_file[:-4]))

    os.system("curl "
              "-F _ggi_project=PPHWeb2 "
              "-F _ggi_origin=query "
              "-F _ggi_target_pipeline=1 "
              "-F MODELNAME=HumVar "  # HUMDIV or HUMVAR
              "-F UCSCDB=hg19 "  # genome
              "-F SNPFUNC=m "  # functional SNP category, c or m
              "-F _ggi_batch_file=@%s "  # filename of variants, format "e.g. chr1:1158631 A/C"
              "-D - http://genetics.bwh.harvard.edu/cgi-bin/ggi/ggi2.cgi "
              "-c %s.humvar.cookies.txt"
              "> %s.humvar.txt" % (input_file, input_file[:-4], input_file[:-4]))

    humdiv_file = "%s.humdiv.txt" % input_file[:-4]
    humvar_file = "%s.humvar.txt" % input_file[:-4]

    return humdiv_file, humvar_file


def get_polyphen_results(hum_file):
    """Gets session ID from input file, gets polyphen results for the session and reads results into a dataframe then
    creates a dictionary (easier to call in template).

    Notes: providing cookies and headers to urllib2 reduces connection errors.

    """
    session_id = None
    opened_file = file(hum_file, "r").read()
    for word in opened_file.split():
        if re.match("polyphenweb2=.{40,};", word):
            session_id = word[13:53]
        else:
            pass
    if session_id is not None:
        polyphen_dict = {}
        cj = cookielib.MozillaCookieJar('%s.cookies.txt' % hum_file[:-4])
        cj.load()
        opener = urllib2.build_opener(urllib2.HTTPCookieProcessor(cj))
        opener.addheaders = [('User-Agent', 'Mozilla/5.0')]
        home = opener.open("http://genetics.bwh.harvard.edu/ggi/pph2/%s/1/pph2-short.txt" % session_id)
        s = home.read()
        with open("hello.txt", "w+") as f:
            f.write(s)
        f.close()
        polyphen_df = pd.read_table("hello.txt", sep='\t', index_col=False)
        polyphen_df = polyphen_df.rename(columns=lambda x: x.strip())
        polyphen_df = polyphen_df.dropna()

        for row_index, row in polyphen_df.iterrows():
            polyphen_dict = {
                'up_acc': row['acc'].strip(),
                'aa1': row['aa1'].strip(),
                'pos': row['pos'],
                'aa2': row['aa2'].strip(),
                'pred': row['prediction'].strip(),
                'prob_score': row['pph2_prob'],
                'sens': row['pph2_TPR'],
                'spec': (100.00 - row['pph2_FPR'])
            }

        return polyphen_dict

    else:
        print "No session ID detected"


def get_polyphen_colormap(prob_score, hum_type):
    """Generates a colour bar similar to that seen in the PolyPhen-2 interface using matplotlib. PNG image saved to
    static folder.

    http://matplotlib.org/examples/api/colorbar_only.html - code here used to help with this.

    """
    fig = pyplot.figure(figsize=(8, 3))
    ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])

    colors = [[0, 255, 0], [2, 255, 0], [4, 255, 0], [6, 255, 0], [8, 255, 0], [10, 255, 0], [12, 255, 0], [14, 255, 0],
              [16, 255, 0], [18, 255, 0], [20, 255, 0], [22, 255, 0], [24, 255, 0], [26, 255, 0], [28, 255, 0],
              [30, 255, 0], [32, 255, 0], [34, 255, 0], [36, 255, 0], [38, 255, 0], [40, 255, 0], [42, 255, 0],
              [44, 255, 0], [46, 255, 0], [48, 255, 0], [50, 255, 0], [52, 255, 0], [54, 255, 0], [56, 255, 0],
              [58, 255, 0], [60, 255, 0], [62, 255, 0], [64, 255, 0], [66, 255, 0], [68, 255, 0], [70, 255, 0],
              [72, 255, 0], [74, 255, 0], [76, 255, 0], [78, 255, 0], [80, 255, 0], [82, 255, 0], [83, 255, 0],
              [85, 255, 0], [87, 255, 0], [89, 255, 0], [91, 255, 0], [93, 255, 0], [95, 255, 0], [97, 255, 0],
              [99, 255, 0], [101, 255, 0], [103, 255, 0], [105, 255, 0], [107, 255, 0], [109, 255, 0], [111, 255, 0],
              [113, 255, 0], [115, 255, 0], [117, 255, 0], [119, 255, 0], [121, 255, 0], [123, 255, 0], [125, 255, 0],
              [127, 255, 0], [129, 255, 0], [131, 255, 0], [133, 255, 0], [135, 255, 0], [137, 255, 0], [139, 255, 0],
              [141, 255, 0], [143, 255, 0], [145, 255, 0], [147, 255, 0], [149, 255, 0], [152, 255, 0], [154, 255, 0],
              [156, 255, 0], [158, 255, 0], [160, 255, 0], [162, 255, 0], [164, 255, 0], [166, 255, 0], [168, 255, 0],
              [170, 255, 0], [172, 255, 0], [174, 255, 0], [176, 255, 0], [178, 255, 0], [180, 255, 0], [182, 255, 0],
              [184, 255, 0], [186, 255, 0], [188, 255, 0], [190, 255, 0], [192, 255, 0], [194, 255, 0], [196, 255, 0],
              [198, 255, 0], [200, 255, 0], [202, 255, 0], [204, 255, 0], [206, 255, 0], [208, 255, 0], [210, 255, 0],
              [212, 255, 0], [214, 255, 0], [216, 255, 0], [218, 255, 0], [220, 255, 0], [222, 255, 0], [224, 255, 0],
              [226, 255, 0], [228, 255, 0], [230, 255, 0], [232, 255, 0], [235, 255, 0], [237, 255, 0], [239, 255, 0],
              [241, 255, 0], [243, 255, 0], [245, 255, 0], [247, 255, 0], [249, 255, 0], [251, 255, 0], [253, 255, 0],
              [255, 255, 0], [255, 255, 0], [255, 253, 0], [255, 251, 0], [255, 249, 0], [255, 247, 0], [255, 245, 0],
              [255, 243, 0], [255, 241, 0], [255, 239, 0], [255, 237, 0], [255, 235, 0], [255, 232, 0], [255, 230, 0],
              [255, 228, 0], [255, 226, 0], [255, 224, 0], [255, 222, 0], [255, 220, 0], [255, 218, 0], [255, 216, 0],
              [255, 214, 0], [255, 212, 0], [255, 210, 0], [255, 208, 0], [255, 206, 0], [255, 204, 0], [255, 202, 0],
              [255, 200, 0], [255, 198, 0], [255, 196, 0], [255, 194, 0], [255, 192, 0], [255, 190, 0], [255, 188, 0],
              [255, 186, 0], [255, 184, 0], [255, 182, 0], [255, 180, 0], [255, 178, 0], [255, 176, 0], [255, 174, 0],
              [255, 172, 0], [255, 170, 0], [255, 168, 0], [255, 166, 0], [255, 164, 0], [255, 162, 0], [255, 160, 0],
              [255, 158, 0], [255, 156, 0], [255, 154, 0], [255, 152, 0], [255, 149, 0], [255, 147, 0], [255, 145, 0],
              [255, 143, 0], [255, 141, 0], [255, 139, 0], [255, 137, 0], [255, 135, 0], [255, 133, 0], [255, 131, 0],
              [255, 129, 0], [255, 127, 0], [255, 125, 0], [255, 123, 0], [255, 121, 0], [255, 119, 0], [255, 117, 0],
              [255, 115, 0], [255, 113, 0], [255, 111, 0], [255, 109, 0], [255, 107, 0], [255, 105, 0], [255, 103, 0],
              [255, 101, 0], [255, 99, 0], [255, 97, 0], [255, 95, 0], [255, 93, 0], [255, 91, 0], [255, 89, 0],
              [255, 87, 0], [255, 85, 0], [255, 83, 0], [255, 82, 0], [255, 80, 0], [255, 78, 0], [255, 76, 0],
              [255, 74, 0], [255, 72, 0], [255, 70, 0], [255, 68, 0], [255, 66, 0], [255, 64, 0], [255, 62, 0],
              [255, 60, 0], [255, 58, 0], [255, 56, 0], [255, 54, 0], [255, 52, 0], [255, 50, 0], [255, 48, 0],
              [255, 46, 0], [255, 44, 0], [255, 42, 0], [255, 40, 0], [255, 38, 0], [255, 36, 0], [255, 34, 0],
              [255, 32, 0], [255, 30, 0], [255, 28, 0], [255, 26, 0], [255, 24, 0], [255, 22, 0], [255, 20, 0],
              [255, 18, 0], [255, 16, 0], [255, 14, 0], [255, 12, 0], [255, 10, 0], [255, 8, 0], [255, 6, 0],
              [255, 4, 0], [255, 2, 0], [255, 0, 0]]

    colour_codes = [[float(j) / 255 for j in i] for i in colors]
    my_cmap = ListedColormap(colour_codes)
    norm = mpl.colors.Normalize(vmin=0.00, vmax=1.00)  # normalise colors to data

    mpl.colorbar.ColorbarBase(ax1, cmap=my_cmap, norm=norm, orientation='horizontal')
    pyplot.plot((prob_score, prob_score), (0.0, 1.0), 'k', linewidth=3)

    fig.savefig("%s_%s.png" % (hum_type, prob_score))
    os.system("mv %s_%s.png %s/static/aml/" % (script_dir, hum_type, prob_score))


