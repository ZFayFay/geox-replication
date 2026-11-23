# I ran 6,300 fake geo experiments on real sales data. Google's placebo test found a lift that wasn't there 1 time in 5.

## Rebuilding Meridian GeoX on ten years of Iowa county sales, with my predictions written down first.

*Part 1 of Big Tech measurement, rebuilt.*

A geo experiment is the measurement tool you reach for when you can't track individual people. You turn ads on in some regions, leave others alone, and compare. I ran them on the agency side, and the hard part was never the arithmetic. It was deciding how much to trust the answer, because regions are few, noisy and all different, and a single bad week in one big city can look like a campaign effect.

In September, Google released [Meridian GeoX](https://research.google/pubs/pub1090193/), an open-source framework for designing and analyzing these tests, with [a paper](https://storage.googleapis.com/gweb-research2023-media/pubtools/1043733.pdf) by Steven Ye, Fabian Sinn, Yunxiao Li and Yang Jiao. It does three things. It picks test and control regions by stratified sampling, grouping regions with similar sales patterns and drawing from each group. It estimates the counterfactual with time-based regression (TBR), synthetic control or synthetic difference-in-differences. And it judges significance with what the authors call design-aware placebo inference: it builds fake treatment groups out of the control regions, using the same rules that built the real one, and asks whether the real result stands out from the fakes. The paper reports that this keeps false positives at the nominal level while cutting the minimum detectable effect by about 30% compared with random placebos. GeoX results are also meant to become priors in Meridian, the MMM I used at Popeyes, which is why I cared about getting this right.

## How I set it up

You can't check a measurement tool on data where you don't know the answer. So I used real sales with fake experiments on top. Iowa publishes every liquor sale from the state wholesaler to stores, back to 2012, and I summed it to weekly sales by county. That gave me 63 studies, each a 50-week window with 26 weeks to train, 8 to validate, 8 to evaluate the design and 8 to run the test, using between 36 and 95 counties depending on which ones had complete data.

Within each study I ran 100 A/A tests, where nothing happens in the test weeks, and I injected known lifts into the treated counties to see whether they were found. Before running anything I wrote down eight hypotheses with pass or fail rules, and committed them to git so the date shows they came first. I rebuilt every step myself from the paper, and then, because the main result surprised me, I ran Google's own library on the same data.

## What held up

The estimators did their job. When I injected a 5% lift, the average estimate was off by only 0.1 to 0.2 percentage points. TBR was also the most efficient of the three estimators, as the paper says. Its median detectable effect was 1.7%, against 2.5% for SDID and 2.7% for synthetic control, and it was the best choice in 43 of 63 studies.

Stratified sampling helped, but less than I expected. It beat random sampling in 38 of 63 studies with a 12% smaller median detectable effect, which missed my pre-registered bar by a hair (p = 0.065). The paper also says synthetic control does better when data is noisy. In my studies the relationship went the other way, though not significantly.

![Design results](https://raw.githubusercontent.com/ZFayFay/geox-replication/main/figures/fig2_design.png)

## What didn't: the uncertainty

The pre-registered test I cared most about was the false positive rate. At a 10% significance level, an honest test flags about 10% of A/A experiments. Design-aware placebo inference, as I implemented it, flagged 22%. A plain t-test that assumes weeks are independent flagged 26%, which is the failure the paper warns about. The surprise was that the simple version Google says it improves on, random placebos, came closest at 12%.

I assumed I had a bug. My first theory was that I was choosing the placebos more strictly than the designs, so I matched the selection exactly. It made no difference (23%). My second was that the placebo groups were smaller than the real ones, so I matched their size. Also no difference. Then I looked at the spread directly. Within a study, the real A/A results varied with a standard deviation of 3.4 on the scale the test uses. The design-aware placebos varied by 2.6, and the random ones by 4.4. The fake groups built with the same rules came out more similar to each other than real designs are, so the yardstick was too short, and ordinary noise looked like an effect.

So I installed Google's `meridian-geox` package and ran it with its defaults on the same Iowa data at daily grain, which is what it expects. Over 60 windows, it flagged 42% of A/A tests at the 10% level. Its estimates for a 5% lift averaged 5.3%, so again the estimate was fine, but its 90% interval contained the true lift only 58% of the time.

![False positives](https://raw.githubusercontent.com/ZFayFay/geox-replication/main/figures/fig1_false_positives.png)

This also changes how to read the claimed power gain. In my runs, design-aware placebos did cut the detectable effect, by 42%, more than the paper's 31%. But a test looks more powerful when its error bars are too narrow, and that is what happened here.

## The design that promised too much

My last test was the one I expected to break. GeoX picks the design with the smallest detectable effect out of many candidates, and the best of many is usually lucky. So I injected a lift exactly equal to the detectable effect each chosen design promised, which should be found 80% of the time. It was found 33% of the time. The chosen designs promised a median of 1.6% and, measured in the actual test weeks, delivered 3.7%. Typical designs, not picked as the best, did better at 47%, but still well short of 80%.

![MDE promise](https://raw.githubusercontent.com/ZFayFay/geox-replication/main/figures/fig3_mde_promise.png)

## What this does and doesn't show

This is one dataset. Iowa's numbers are orders from stores, which are lumpy, especially by day, and counties aren't the media markets advertisers use. Google tested on 192 advertisers' data and reported clean results. My rebuild also simplifies parts of their design search. So I am not saying GeoX is broken. I am saying that on a real, public dataset where anyone can check my work, both my implementation and theirs were far more confident than they should have been, and the estimates were fine while the intervals weren't.

## What I take from it

Before trusting any geo test tool, including this one, I would backtest it on my own history. Take two years of sales, pretend a campaign ran in a dozen past periods when nothing happened, and count how often the tool says something did. That costs an afternoon, and it tells you how wide your error bars really are. It matters even more if the result becomes a prior in an MMM, because an interval that's too tight there quietly drags every channel's ROI toward the test. And when a tool picks the best of many designs, I'd plan the budget around a detectable effect at least twice as big as the one it promises.

---

*Notebook, data and code: [github.com/ZFayFay/geox-replication](https://github.com/ZFayFay/geox-replication). The pre-registration is the first commit, and the library check is in `libcheck/`. GeoX belongs to Steven Ye, Fabian Sinn, Yunxiao Li, Yang Jiao and the Google Meridian team. This is an independent rebuild of their [paper](https://research.google/pubs/pub1090193/) using their [open-source code](https://github.com/google/meridian-geox).*
