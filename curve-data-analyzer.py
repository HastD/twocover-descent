import json

with open("data/labels.json", "r") as f:
    labels = json.load(f)

curves = []
logs = []
for label in labels:
    with open("../twocover-results/json/curve-{}.json".format(label), "r") as f:
        curves.append(json.load(f))
    with open("../twocover-results/logs/curve-{}.log".format(label), "r") as f:
        logs.append(f.read())

# Curves where the method succeeded
success = [C for C in curves if C["verified"]]
print("succeeded: {}".format(len(success)))

remaining = [C for C in curves if C not in success]

# Curves where a failure of the Hasse principle is an obstruction
hasse = [C for C in remaining if C["hasse_principle"] is False]
print("Hasse principle failure: {}".format(len(hasse)))

remaining = [C for C in remaining if C not in hasse]

# Curves where a MW rank being too high is an obstruction
high_rank = [C for C in remaining if any(
    not twist["verified"]
    and any(D["rank"] is not None and D["rank"] >= len(D["g"]) - 1 for D in twist["g1"])
    for twist in C["twists"])]
print("high rank: {}".format(len(high_rank)))

# Curves where Magma failing to provably compute a MW group is an obstruction
no_MW = [C for C in remaining if any(
    not twist["verified"]
    and any(D["MW_proven"] is False for D in twist["g1"])
    for twist in C["twists"])]
print("can't compute MW: {}".format(len(no_MW)))

remaining = [C for C in remaining if C not in high_rank and C not in no_MW]

bins = 10
disc_bin = lambda C, i: (i*1e5 < C["disc"] <= (i + bins/10)*1e5)

success_runtimes_vs_disc_median = [median([C["runtime"] for C in success if disc_bin(C, i)]) for i in range(bins)]
success_runtimes_vs_disc_mean = [mean([C["runtime"] for C in success if disc_bin(C, i)]) for i in range(bins)]

success_rates_vs_disc = [len([C for C in success if disc_bin(C, i)])/len([C for C in curves if disc_bin(C, i)]) for i in range(bins)]
hasse_rates_vs_disc = [len([C for C in hasse if disc_bin(C, i)])/len([C for C in curves if disc_bin(C, i)]) for i in range(bins)]
high_rates_vs_disc = [len([C for C in high_rank if disc_bin(C, i)])/len([C for C in curves if disc_bin(C, i)]) for i in range(bins)]
noMW_rates_vs_disc = [len([C for C in no_MW if disc_bin(C, i)])/len([C for C in curves if disc_bin(C, i)]) for i in range(bins)]
timeout_rates_vs_disc = [len([C for C in remaining if disc_bin(C, i)])/len([C for C in curves if disc_bin(C, i)]) for i in range(bins)]

plot_disc_rates = lambda rates: [((i + bins/20)*1e5, rates[i]) for i in range(bins)]

list_plot(plot_disc_rates(success_rates_vs_disc), plotjoined=True, xmin=0, xmax=1e6, ymin=0, color="blue") \
+ list_plot(plot_disc_rates(hasse_rates_vs_disc), plotjoined=True, xmin=0, xmax=1e6, ymin=0, color="red")
# + list_plot(plot_disc_rates(high_rates_vs_disc), plotjoined=True, color="purple") \
# + list_plot(plot_disc_rates(noMW_rates_vs_disc), plotjoined=True, color="red") \
# + list_plot(plot_disc_rates(timeout_rates_vs_disc), plotjoined=True, color="blue")

