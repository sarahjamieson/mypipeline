import os

script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(script_dir, os.pardir))


def generate_ngs_vs_frag(frag, ngs, itd):
    with open("all_samples.R", "w+") as r_file:
        r_file.write("library(ggplot2)\npng(\"%s/static/aml/all_samples.png\")\nFrag_Analysis <- c(" % script_dir)
        for f in frag[:-1]:
            r_file.write("%s, " % f)
        r_file.write("%s)\nNGS <- c(" % frag[-1])
        for n in ngs[:-1]:
            r_file.write("%s, " % n)
        r_file.write("%s)\nITD_size <- c(" % ngs[-1])
        for s in itd[:-1]:
            r_file.write("%s, " % s)
        r_file.write("%s)\ndf <- data.frame(Frag_Analysis, NGS, ITD_size)\n"
                     "p <- ggplot(data = df, aes(x = Frag_Analysis, NGS))\n"
                     "p + geom_point(aes(size=ITD_size)) + xlab(\"Fragment Analysis AR\") + "
                     "ylab(\"NGS AR\")\ndev.off()\n" % itd[-1])
    r_file.close()
    os.system("Rscript %s/all_samples.R" % parent_dir)


def generate_ngs_vs_frag_exclude_outliers(frag, ngs, itd):
    with open("all_samples_no_outliers.R", "w+") as r_file:
        r_file.write(
            "library(ggplot2)\npng(\"%s/static/aml/all_samples_no_outliers.png\")\nFrag_Analysis <- c(" % script_dir)
        for f in frag[:-1]:
            if f > 5.00:
                # get_index and remove from ngs and itd at that index
                f_index = frag.index(f)
                ngs.pop(f_index)
                itd.pop(f_index)
            else:
                r_file.write("%s, " % f)
        r_file.write("%s)\nNGS <- c(" % frag[-1])
        for n in ngs[:-1]:
            r_file.write("%s, " % n)
        r_file.write("%s)\nITD_size <- c(" % ngs[-1])
        for s in itd[:-1]:
            r_file.write("%s, " % s)
        r_file.write("%s)\ndf <- data.frame(Frag_Analysis, NGS, ITD_size)\n"
                     "p <- ggplot(data = df, aes(x = Frag_Analysis, NGS))\n"
                     "p + geom_point(aes(size=ITD_size)) + xlab(\"Fragment Analysis AR\") + "
                     "ylab(\"NGS AR\")\ndev.off()\n" % itd[-1])
    r_file.close()
    os.system("Rscript %s/all_samples_no_outliers.R" % parent_dir)

