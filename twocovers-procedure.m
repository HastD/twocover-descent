/*
Provably computes rational points on a curve using two-cover descent and elliptic Chabauty.
*/

// Genus 2 curve 6443.a.6443.1
R<x> := PolynomialRing(Rationals());
C := HyperellipticCurve(R![0, 1, 1, -2, -1, 1], R![1]); C;

load "twocovers.m";

SetVerbose("MWSha", 1);
SetVerbose("MordellWeilGroup", 1);
SetVerbose("Selmer", 1);
SetVerbose("TwoDescent", 1);
SetVerbose("CasselsTate", 1);
SetVerbose("EllChab", 1);
//SetVerbose("ClassGroup", 3);
SetVerbose("LocSol", 1);

SetLogFile("twocover-descent-log.txt");

descent_procedure(C : SearchBound := 10000, AssumeGRH := true, OutputFile := "Output/twocovers-6443.a.6443.1-GRH-results.txt");

