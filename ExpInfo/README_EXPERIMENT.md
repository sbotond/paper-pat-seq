## RNA isolation

Wt yeast (BY4741, which is derived from S288c) and Δccr4/pan2 (Beilharz et 
al., RNA. 2007 Jul;13(7):982-97) were single colony streaked on YPAD-agar 
plates, and then 2 individual colonies each were inoculated in 2ml YPAD medium 
and grown o.n. at 30°C at 180rpm. OD was measured and 7x107 cells were 
pelleted, resuspended in 2ml Y1 medium (2μl beta-ME, 150μl zymolase) and 
incubated for 30min at 30°C in a waterbath. After centrifugation, the pellet 
was lysed in buffer RLT and processed according to the yeast protocol of the 
Qiagen RNeasy kit. Yield and quality of totalRNA was measured with Nanodrop and 
confirmed by Agilent RNA Nano chip (RIN 10). A DNase digestion was performed on 
20μg totalRNA using Turbo DNase (Life Technologies), followed by an RNeasy 
column cleanup step and elution into 30μl nuclease free water (Ambion).

## PolyA selection

2 rounds of polyA selection were performed using the Truseq RNA Sample 
Preparation kit (Illumina) with several modifications. In more detail: 
50μl of DNase digested totalRNA was incubated with 50μl RNA Purification 
beads (oligo-dT) for 5min at 65°C, cooled to 4°C, and then incubated for 5min 
at RT in a DNA Engine tetrade 2 (Biorad) thermocycler. Afterwards, beads were 
separated with a magnet and washed with 200μl Bead Wash buffer, and then 
eluted in 50μl elution buffer for 2min at 80°C. On RT, 50μl of binding 
buffer were added, placed on a magnet and washed with 200μl Bead Wash buffer, 
and again eluted in 50μl elution buffer for 2min at 80°C. The supernatant 
containing the mRNA was transferred into a 1.5ml Low Binding tube (Eppendorf) 
and quantified with the Qubit RNA assay.

## Spike-in control

20μl of spike in in vitro transcripts with a known A42 stretch (Pelechano V. 
et al., Nature. 2013 May 2;497(7447):127-31) were added. A cleanup step was 
performed using the RNeasy MinElute Cleanup Kit (Qiagen) according to protocol 
eluting with 14μl water.

## G-tailing and fragmentation

The USB Poly(A) Tail-Length Assay kit (Affymetrix) was used with several 
modifications to add a short G stretch to the 3’end of the mRNA. In detail, 
14μl mRNA (with spike-in) were incubated with 4μl of 5xTail buffer mix and 
2μl 10x Tail Enzyme mix at 37°C for 1h and stopped by addition of 2μl 10x 
Tail stop solution. After an RNeasy Minelute cleanup step and elution into 
18μl of water, chemical fragmentation using the NEBnext Magnesium RNA 
Fragmentation kit (NEB) was performed in 20μl for 3.5min at 94°C. Afterwards, 
the reaction was immediately put on ice, and 2μl of RNA Fragmentation Stop 
Solution was added. Another RNeasy Minelute cleanup step was performed, eluting 
into 15μl of water.

## cDNA synthesis

1ul TTTTVN (16.7μM) and 1ul CCCCCCTT (50μM) custom primers (both Sigma) were 
added to the 15μl fragmented RNA, primed by incubation at 65°C for 5min and 
put on ice. 1st strand cDNA synthesis was performed using the 1st strand 
mastermix (Illumina Truseq RNA kit) and SuperscriptII enzyme (Life 
Technologies) by incubating at 25°C for 10min, 42°C for 50min at 70°C for 
15min, and cooled to 4°C in a MJ-Mini thermocycler (Biorad). 2nd strand 
synthesis was performed by addition of 25μl 2nd strand master mix and 
incubation for 1h at 16°C.

## NGS library preparation

An AMPure XP (Beckman Coulter) cleanup step was performed by adding 90μl 
AMPure XP beads, vortexing and incubating for 15min at RT. After separation on 
a magnet, the supernatant was removed and the beads were washed twice with 
200μl 80% ethanol. After air drying for 10min at RT, cDNA was eluted in 50μl 
Resuspension buffer and stored at -20°C.
End repair was performed by adding 10μl Resuspension buffer and 40μl 
Endrepair mix and incubating at 30°C for 30min at 750rpm in a thermomixer 
(Eppendorf).
After another AMPure XP cleanup step using 160μl beads and eluting into 
17.5μl Resuspension buffer, the A-overhang addition step was performed. 
12.5μl A-tailing mix were added and incubated at 37°C for 30min at 750rpm, 
followed by 5min at 70°C and put on ice. Adapter ligation used 2.5μl of a 
1:40 dilution of Adapter oligo mix (Agilent Sureselect v4 kit), 2.5μl 
Resuspension solution and 2.5μl DNA ligase mix, incubated at 30°C for 10min 
in a water bath and stopped by addition of 5μl Stop ligase mix. Adapter 
ligated cDNA library was cleaned up by 2 rounds of AMPure XP using first 42μl 
of AMPure XP beads and eluting in 50μl Resuspension buffer and then using 
50μl AMPure XP beads and eluting in 40μl Resuspension buffer.
Next, a size selection step was performed to create 4 technical replicates 
between 300-325bp, 325-350bp, 400-425bp and 425-450bp. For this, the purified 
adapter ligated cDNA was loaded on a 1.5% Agarose gel and separated for 2h at 
120V. Using SybrSafe (Invitrogen) as dye and a Dark Reader transilluminator 
(Clare Chemical Research), narrow bands were cut using 100bp marker as guidance 
and gel extracted using the Nucleospin Gel and PCR Cleanup kit (Macherey-Nagel) 
and eluted in 22μl NE buffer.
Two sequential rounds of PCR were performed to introduce barcodes and amplify 
the cDNA library. First, the Sureselect primer and Sureselect ILM Indexing 
Pre-Capture PCR reverse primer were used with the Phusion HF 2x mix in 50μl in 
a MJ Mini thermocycler as follows: 2min 98°C, followed by 5 cycles of 30sex 
98°C, 30sec 65°C, 1min at 72°C, and a final 72°C step for 10min. Next, an 
AMPure XP cleanup step was performed using 90μl AMPure XP beads and eluting 
into 23μl water. Afterwards, a 2nd PCR was performed using the Sureselect ILM 
Indexing Post Capture Forward PCR primer and a different PCR Primer Index 
primer/sample together with the Phusion HF 2x mix in 50μl in a MJ Mini 
thermocycler as follows: 2min 98°C, followed by 13-15 cycles of 30sex 98°C, 
30sec 57°C, 1min at 72°C, and a final 72°C step for 10min. After a final 
AMPure XP cleanup step using 90μl beads, the final library was eluted in 20μl 
water, quantified with the Qubit HS dsDNA kit and pooled equimolarly.

## NGS Sequencing

All 16 pooled barcoded libraries were sequenced on 1 Hiseq2000 lane in 2x100bp 
mode.

