import QtQuick 2.9
import QtQuick.Controls 2.2
import QtQuick.Layouts 1.3
import QtQuick.Dialogs 1.2
import io.qt.examples.backend 1.0
import QtQuick.Controls.Styles 1.4
import QtQuick.Window 2.2
import QtCharts 2.2

import FileIO 1.0

ApplicationWindow {
    id: root
    visible: true
    width: 1440
    height: 900
    title: qsTr("MCMD :: Monte Carlo & Molecular Dynamics")

    // custom c++ classes
    BackEnd {
        id: backend
        outputLineNumber: 0
    }
    FileIO {
        id: runlogFile
        property int colOffset: 2; // to read data excluding the line # output
        //homeDir: "/Users/douglasfranz" // determined dynamically now
        //source: homeDir+"/mcmd/testzone/runlog.log" // determined dynamically now
        onError: console.log(msg);
    }

    // main visual stuff begins
    SwipeView {
        id: swipeView
        anchors.fill: parent
        currentIndex: tabBar.currentIndex
        onCurrentIndexChanged: {
            //myText.text = myFile.read();
        }

        Page { // 1 :: Input stuff
            Rectangle {
                id: leftcol
                color: "red"
                border.color: "black"
                border.width: 3
                height: parent.height
                width: parent.width/2
                anchors.margins: 5

               Rectangle {
                    width: parent.width
                    height: 50
                    color: "#cdcdcd"
                    Text {
                        font.pixelSize: 30
                        text: "MCMD input parameters"
                    }
                }
/*
                    InputLine {
                        id: firstInpLine
                        y: 50
                        name: "Job name"
                        defaultIn: "test job"
                    }
                    InputLine {
                        y: firstInpLine.y + firstInpLine.height
                        name: "Mode"
                        defaultIn: "mc"
                    }
                    InputLine {
                        y: firstInpLine.y + 2*firstInpLine.height
                        name: "Input atoms"
                        defaultIn: "/Users/douglasfranz/mcmd/testzone/input.pdb"
                    }
                    InputLine {
                        y: firstInpLine.y + 3*firstInpLine.height
                        name: "Potential"
                        defaultIn: "ljes"
                    }
                    InputLine {
                        y: firstInpLine.y + 4*firstInpLine.height
                        name: "Sorbates"
                        defaultIn: "h2_bss"
                    }
                    InputLine {
                        y: firstInpLine.y + 5*firstInpLine.height
                        name: "Fugacity"
                        defaultIn: "h2"
                    }
                    InputLine {
                        y: firstInpLine.y + 6*firstInpLine.height
                        name: "Basis a"
                        defaultIn: "25.669"
                    }
                    InputLine {
                        y: firstInpLine.y + 7*firstInpLine.height
                        name: "Basis b"
                        defaultIn: "25.669"
                    }
                    InputLine {
                        y: firstInpLine.y + 8*firstInpLine.height
                        name: "Basis c"
                        defaultIn: "25.669"
                    }
                    InputLine {
                        y: firstInpLine.y + 9*firstInpLine.height
                        name: "Basis α"
                        defaultIn: "90.0"
                    }
                    InputLine {
                        y: firstInpLine.y + 10*firstInpLine.height
                        name: "Basis β"
                        defaultIn: "90.0"
                    }
                    InputLine {
                        y: firstInpLine.y + 11*firstInpLine.height
                        name: "Basis γ"
                        defaultIn: "90.0"
                    }
                    InputLine {
                        y: firstInpLine.y + 12*firstInpLine.height
                        name: "Ensemble"
                        defaultIn: "uvt"
                    }
                    InputLine {
                        y: firstInpLine.y + 13*firstInpLine.height
                        name: "Corrtime"
                        defaultIn: "1000"
                    }
                    InputLine {
                        y: firstInpLine.y + 14*firstInpLine.height
                        name: "Final step"
                        defaultIn: "10000000"
                    }
                    InputLine {
                        y: firstInpLine.y + 15*firstInpLine.height
                        name: "Temperature"
                        defaultIn: "77.0"
                    }
                    InputLine {
                        y: firstInpLine.y + 16*firstInpLine.height
                        name: "Pressure"
                        defaultIn: "1.0"
                    }
                    InputLine {
                        y: firstInpLine.y + 17*firstInpLine.height
                        name: "Insert factor"
                        defaultIn: "0.667"
                    }
                    InputLine {
                        y: firstInpLine.y + 18*firstInpLine.height
                        name: "Displace factor"
                        defaultIn: "2.5"
                    }
                    InputLine {
                        y: firstInpLine.y + 19*firstInpLine.height
                        name: "Angle rotation factor"
                        defaultIn: "360.0"
                    }
                    InputLine {
                        y: firstInpLine.y + 20*firstInpLine.height
                        name: "Basis γ"
                        defaultIn: "90.0"
                    }
                    InputLine {
                        y: firstInpLine.y + 21*firstInpLine.height
                        name: "Feynman-Hibbs corrections"
                        defaultIn: "on"
                    }
                    InputLine {
                        y: firstInpLine.y + 22*firstInpLine.height
                        name: "F-H order"
                        defaultIn: "4"
                    }
                    InputLine {
                        y: firstInpLine.y + 23*firstInpLine.height
                        name: "write_lammps"
                        defaultIn: "on"
                    }
                    InputLine {
                        y: firstInpLine.y + 24*firstInpLine.height
                        name: "Auto-reject option"
                        defaultIn: "on"
                    }
                    InputLine {
                        y: firstInpLine.y + 25*firstInpLine.height
                        name: "Auto-reject r"
                        defaultIn: "1.78"
                    }
*/
            }
            Rectangle {
                id: midcol
                anchors.left: leftcol.right
                border.color: "black"
                border.width: 3
                color: "blue"
                height: parent.height
                width: parent.width/4
            }
            Rectangle {
                id: rightcol
                anchors.left: midcol.right
                color: "orange"
                border.color: "black"
                border.width: 3
                height: parent.height
                width: parent.width/4
            }
            Timer {
                id: timer
                //property BackEnd backend: BackEnd {}
                interval: 300 // ms
                repeat: true
                running: root.visible
                onTriggered: {

                    var newText = runlogFile.read();

                    //scroller.position = (outputText.contentHeight - outputPage.height)/outputText.contentHeight;

                    if (newText !== "") { // make sure we have new data to work with.
                        var steps = new Array();
                        var progresss = new Array();
                        // MC...
                        var PolEns = new Array();
                        var ESEns = new Array();
                        var RDEns = new Array();
                        var UEns = new Array();
                        var Qsts = new Array();
                        var Ns = new Array();
                        var Navgs = new Array();
                        var Inserts = new Array();
                        var Removes = new Array();
                        var Displaces = new Array();
                        var VolumeChanges = new Array();

                        // MD...
                        var KEs = new Array();
                        var PEs = new Array();
                        var TEs = new Array();
                        var Ds = new Array();
                        var ETemps = new Array();
                        var ITemps = new Array();
                        var Press = new Array();
                        var Vels = new Array();


                        // loop each NEW line...
                        var allLines = newText.split("\n"); //outputText.text.split("\n");
                        var i=0;
                        while (allLines.length > i) {
                            var thisLine = allLines[i].split(/(\s)/).filter( function(e) { return e.trim().length > 0; } );;
                            if (allLines[i].indexOf("Step") !== -1) {
                                var step = thisLine[1+runlogFile.colOffset];
                                steps.push(step);
                                var prg = thisLine[6+runlogFile.colOffset].replace("%","");
                                progresss.push(parseFloat(prg));
                            }
                            // MD...
                            else if (allLines[i].indexOf("KE") !== -1) {
                                var KE = thisLine[2+runlogFile.colOffset];
                                KEs.push(KE);
                            }
                            else if (allLines[i].indexOf("PE") !== -1) {
                                var PE = thisLine[2+runlogFile.colOffset];
                                PEs.push(PE);
                            }
                            else if (allLines[i].indexOf("Total E") !== -1) {
                                var TE = thisLine[3+runlogFile.colOffset];
                                TEs.push(TE);
                            } else if (allLines[i].indexOf("Diffusion c") !== -1) {
                                var Diff = thisLine[5+runlogFile.colOffset];
                                Ds.push(Diff);
                            } else if (allLines[i].indexOf("Emergent T") !== -1) {
                                var ETemp = thisLine[3+runlogFile.colOffset];
                                ETemps.push(ETemp);
                            } else if (allLines[i].indexOf("Instantaneous T") !== -1) {
                                var ITemp = thisLine[3+runlogFile.colOffset];
                                ITemps.push(ITemp);
                            } else if (allLines[i].indexOf("Pressure") !== -1) {
                                var Pres = thisLine[3+runlogFile.colOffset];
                                Press.push(Pres);
                            } else if (allLines[i].indexOf("Average v =") !== -1) {
                                var Vel = thisLine[3+runlogFile.colOffset];
                                Vels.push(Vel);
                            // MC ...
                            } else if (allLines[i].indexOf("Polar avg") !== -1) {
                                var PolEn = thisLine[3+runlogFile.colOffset];
                                PolEns.push(PolEn);
                            } else if (allLines[i].indexOf("ES avg =") !== -1) {
                                var ESEn = thisLine[3+runlogFile.colOffset];
                                ESEns.push(ESEn);
                            } else if (allLines[i].indexOf("RD avg =") !== -1) {
                                var RDEn = thisLine[3+runlogFile.colOffset];
                                RDEns.push(RDEn);
                            } else if (allLines[i].indexOf("Total potential avg") !== -1) {
                                var UEn = thisLine[4+runlogFile.colOffset];
                                UEns.push(UEn);
                            } else if (allLines[i].indexOf("Qst = ") !== -1) {
                                var Qst = thisLine[2+runlogFile.colOffset];
                                Qsts.push(Qst);
                            } else if (allLines[i].indexOf("N_movables =") !== -1 && allLines[i].indexOf("mg/g") !== -1) {
                                var Navg = thisLine[2+runlogFile.colOffset];
                                Navgs.push(Navg);
                            } else if (allLines[i].indexOf("N_movables =") !== -1 && allLines[i].indexOf("mg/g") === -1) {
                                var N = thisLine[5+runlogFile.colOffset];
                                Ns.push(N);
                            } else if (allLines[i].indexOf("Total accepts:") !== -1) {
                                var ins = thisLine[4+runlogFile.colOffset].replace("%","");
                                Inserts.push(ins);
                                var rem = thisLine[7+runlogFile.colOffset].replace("%","");
                                Removes.push(rem);
                                var dis = thisLine[10+runlogFile.colOffset].replace("%","");
                                Displaces.push(dis);
                                var vc = thisLine[13+runlogFile.colOffset].replace("%","");
                                VolumeChanges.push(vc);
                            }
                            i++;
                        }

                        var laststep = steps[steps.length-1];
                        progressdisplayer.update(progresss[progresss.length-1]);
                        mdprogressdisplayer.update(progresss[progresss.length-1]);
                        energychart.updateKE(laststep, KEs[KEs.length -1]);
                        energychart.updatePE(laststep, PEs[PEs.length -1]);
                        energychart.updateTE(laststep, TEs[TEs.length -1]);
                        diffusionchart.updateDiff(laststep, Ds[Ds.length -1]);
                        temperaturechart.updateETemp(laststep, ETemps[ETemps.length -1]);
                        temperaturechart.updateITemp(laststep, ITemps[ITemps.length -1]);
                        pressurechart.updatePres(laststep, Press[Press.length -1]);
                        velocitychart.updateVel(laststep, Vels[Vels.length -1]);
                        mcenergychart.updatePolEn(laststep, PolEns[PolEns.length -1]);
                        mcenergychart.updateESEn(laststep, ESEns[ESEns.length -1]);
                        mcenergychart.updateRDEn(laststep, RDEns[RDEns.length -1]);
                        mcenergychart.updateTotalU(laststep, UEns[UEns.length -1]);
                        qstchart.updateQst(laststep, Qsts[Qsts.length -1]);
                        nchart.updateNAvg(laststep, Navgs[Navgs.length -1]);
                        nchart.updateN(laststep, Ns[Ns.length -1]);
                        mcstatschart.updateMCStats(Inserts[Inserts.length-1],
                            Removes[Removes.length-1], Displaces[Displaces.length-1],
                            VolumeChanges[VolumeChanges.length-1]);
                    } // end if new data is detected
                }
            } // end timer

        }
        Page { // 2 :: MC graphs
            id: mcgraphspage
            Text {
                id: mctoptitle
                text: "Live data :: progress: "
                height: 25
                width: 150

            }
            ProgressBar {
                id: progressdisplayer
                height: 25
                anchors.left: mctoptitle.right
                value: 0
                function update(newprog) {
                    this.value = newprog;
                }
            }

            ChartView {
                id: mcenergychart
                theme: ChartView.ChartThemeDark
                title: "Potential Energy"

                //anchors.fill: parent
                anchors.left: parent.left
                anchors.top: mctoptitle.bottom
                height: (root.height - mctoptitle.height)/2
                width: root.width/3
                antialiasing: true
                function updatePolEn(x, y) {
                    pol_line.append(x,y);
                    if (x > pol_line.axisX.max) {
                        pol_line.axisX.max = x;
                    }
                    else if (x < pol_line.axisX.min) {
                        pol_line.axisX.min = x;
                    }
                    if (y > pol_line.axisY.max) {
                        pol_line.axisY.max = y;
                    }
                    else if (y < pol_line.axisY.min) {
                        pol_line.axisY.min = y;
                    }
                }
                function updateESEn(x, y) {
                    es_line.append(x,y);
                    if (x > es_line.axisX.max) {
                        es_line.axisX.max = x;
                    }
                    else if (x < es_line.axisX.min) {
                        es_line.axisX.min = x;
                    }
                    if (y > es_line.axisY.max) {
                        es_line.axisY.max = y;
                    }
                    else if (y < es_line.axisY.min) {
                        es_line.axisY.min = y;
                    }
                }
                function updateRDEn(x, y) {
                    rd_line.append(x,y);
                    if (x > rd_line.axisX.max) {
                        rd_line.axisX.max = x;
                    }
                    else if (x < rd_line.axisX.min) {
                        rd_line.axisX.min = x;
                    }
                    if (y > rd_line.axisY.max) {
                        rd_line.axisY.max = y;
                    }
                    else if (y < rd_line.axisY.min) {
                        rd_line.axisY.min = y;
                    }
                }
                function updateTotalU(x, y) {
                    totalu_line.append(x,y);
                    if (x > totalu_line.axisX.max) {
                        totalu_line.axisX.max = x;
                    }
                    else if (x < totalu_line.axisX.min) {
                        totalu_line.axisX.min = x;
                    }
                    if (y > totalu_line.axisY.max) {
                        totalu_line.axisY.max = y;
                    }
                    else if (y < totalu_line.axisY.min) {
                        totalu_line.axisY.min = y;
                    }
                }

                LineSeries {
                    id: pol_line
                    axisX: ValueAxis {
                        min: 0
                        max: 0
                        tickCount: 5
                        titleText: "Step"
                    }

                    axisY: ValueAxis {
                        min: 0
                        max: 1e-7
                        titleText: "Energy (K)"
                    }
                    name: "Polarization"
                }
                LineSeries {
                    id: es_line
                    name: "Electrostatics"
                }
                LineSeries {
                    id: rd_line
                    name: "Repulsion/dispersion"
                }
                LineSeries {
                    id: totalu_line
                    name: "Total U"
                }
            }
            ChartView {
                id: qstchart
                theme: ChartView.ChartThemeDark
                title: "Isosteric Heat"

                //anchors.fill: parent
                anchors.left: mcenergychart.right
                anchors.top: mctoptitle.bottom
                height: (root.height - mctoptitle.height)/2
                width: root.width/3
                antialiasing: true
                function updateQst(x, y) {
                    qst_line.append(x,y);
                    if (x > qst_line.axisX.max) {
                        qst_line.axisX.max = x;
                    }
                    else if (x < qst_line.axisX.min) {
                        qst_line.axisX.min = x;
                    }
                    if (y > qst_line.axisY.max) {
                        qst_line.axisY.max = y;
                    }
                    else if (y < qst_line.axisY.min) {
                        qst_line.axisY.min = y;
                    }
                }

                LineSeries {
                    id: qst_line
                    axisX: ValueAxis {
                        min: 0
                        max: 0
                        tickCount: 5
                        titleText: "Step"
                    }

                    axisY: ValueAxis {
                        min: 0
                        max: 1e-7
                        titleText: "Energy (kJ/mol)"
                    }
                    name: "Qst"
                }
            }
            ChartView {
                id: nchart
                theme: ChartView.ChartThemeDark
                title: "Uptake (N)"

                //anchors.fill: parent
                anchors.left: qstchart.right
                anchors.top: mctoptitle.bottom
                height: (root.height - mctoptitle.height)/2
                width: root.width/3
                antialiasing: true
                function updateN(x, y) {
                    n_line.append(x,y);
                    if (x > n_line.axisX.max) {
                        n_line.axisX.max = x;
                    }
                    else if (x < n_line.axisX.min) {
                        n_line.axisX.min = x;
                    }
                    if (y > n_line.axisY.max) {
                        n_line.axisY.max = y;
                    }
                    else if (y < n_line.axisY.min) {
                        n_line.axisY.min = y;
                    }
                }
                function updateNAvg(x, y) {
                    navg_line.append(x,y);
                    if (x > navg_line.axisX.max) {
                        navg_line.axisX.max = x;
                    }
                    else if (x < navg_line.axisX.min) {
                        navg_line.axisX.min = x;
                    }
                    if (y > navg_line.axisY.max) {
                        navg_line.axisY.max = y;
                    }
                    else if (y < navg_line.axisY.min) {
                        navg_line.axisY.min = y;
                    }
                }

                LineSeries {
                    id: n_line
                    axisX: ValueAxis {
                        min: 0
                        max: 0
                        tickCount: 5
                        titleText: "Step"
                    }

                    axisY: ValueAxis {
                        min: 0
                        max: 1e-7
                        titleText: "No. molecules"
                    }
                    name: "N"
                }
                LineSeries {
                    id: navg_line
                    name: "N average"
                }
            }

            ChartView {
                id: mcstatschart
                theme: ChartView.ChartThemeDark
                title: "MC stats"
                anchors.left: parent.left
                anchors.top: mcenergychart.bottom
                height: (root.height - mctoptitle.height)/2
                width: root.width/3
                antialiasing: true
                function updateMCStats(ins, rem, dis, vol) {
                    mcs1.value = ins;
                    mcs2.value = rem;
                    mcs3.value = dis;
                    mcs4.value = vol;

                }

                PieSeries {
                    id: mcstats_pie
                    PieSlice { id: mcs1; color: "green"; label: "Inserts"; value: 94.7 }
                    PieSlice { id: mcs2; color: "red"; label: "Removes"; value: 5.1 }
                    PieSlice { id: mcs3; color: "grey"; label: "Displaces"; value: 0.1 }
                    PieSlice { id: mcs4; color: "blue"; label: "Volume changes"; value: 0.1 }
                }
            }
        }

        Page { // 3 :: MD graphs
            id: mdgraphspage

            Text {
                id: toptitle
                text: "Live data :: progress: "
                height: 25
                width: 150
            }
            ProgressBar {
                id: mdprogressdisplayer
                value: 0
                height: 25

                anchors.left: toptitle.right
                function update(newprog) {
                    this.value = newprog;
                }
            }

            ChartView {
                id: energychart
                theme: ChartView.ChartThemeDark
                title: "Energies"

                //anchors.fill: parent
                anchors.left: parent.left
                anchors.top: toptitle.bottom
                height: (root.height - toptitle.height)/2
                width: root.width/3
                antialiasing: true
                function updateKE(x, y) {
                    ke_line.append(x,y);
                    if (x > ke_line.axisX.max) {
                        ke_line.axisX.max = x;
                    }
                    else if (x < ke_line.axisX.min) {
                        ke_line.axisX.min = x;
                    }
                    if (y > ke_line.axisY.max) {
                        ke_line.axisY.max = y;
                    }
                    else if (y < ke_line.axisY.min) {
                        ke_line.axisY.min = y;
                    }
                }
                function updatePE(x, y) {
                    pe_line.append(x,y);
                    if (x > ke_line.axisX.max) {
                        ke_line.axisX.max = x;
                    }
                    else if (x < ke_line.axisX.min) {
                        ke_line.axisX.min = x;
                    }
                    if (y > ke_line.axisY.max) {
                        ke_line.axisY.max = y;
                    }
                    else if (y < ke_line.axisY.min) {
                        ke_line.axisY.min = y;
                    }
                }
                function updateTE(x, y) {
                    te_line.append(x,y);
                    if (x > ke_line.axisX.max) {
                        ke_line.axisX.max = x;
                    }
                    else if (x < ke_line.axisX.min) {
                        ke_line.axisX.min = x;
                    }
                    if (y > ke_line.axisY.max) {
                        ke_line.axisY.max = y;
                    }
                    else if (y < ke_line.axisY.min) {
                        ke_line.axisY.min = y;
                    }
                }

                LineSeries {
                    id: ke_line
                    axisX: ValueAxis {
                        min: 0
                        max: 10
                        tickCount: 5
                        titleText: "Step"

                    }

                    axisY: ValueAxis {
                        min: -0.5
                        max: 1.5
                        titleText: "Energy (kJ/mol)"
                    }
                    name: "Kinetic Energy"
                }
                LineSeries {
                    id: pe_line
                    name: "Potential Energy"
                }
                LineSeries {
                    id: te_line
                    name: "Total Energy"
                }
            }
            ChartView {
                id: diffusionchart
                theme: ChartView.ChartThemeDark
                title: "Diffusion Coefficient"

                //anchors.fill: parent
                anchors.left: energychart.right
                anchors.top: toptitle.bottom
                height: (root.height - toptitle.height)/2
                width: root.width/3
                antialiasing: true
                function updateDiff(x, y) {
                    diff_line.append(x,y);
                    if (x > diff_line.axisX.max) {
                        diff_line.axisX.max = x;
                    }
                    else if (x < diff_line.axisX.min) {
                        diff_line.axisX.min = x;
                    }
                    if (y > diff_line.axisY.max) {
                        diff_line.axisY.max = y;
                    }
                    else if (y < diff_line.axisY.min) {
                        diff_line.axisY.min = y;
                    }
                }

                LineSeries {
                    id: diff_line
                    axisX: ValueAxis {
                        min: 0
                        max: 0
                        tickCount: 5
                        titleText: "Step"
                    }

                    axisY: ValueAxis {
                        min: 0
                        max: 1e-7
                        titleText: "D (cm^2 / s)"
                    }
                    name: "D"
                }
            }
            ChartView {
                id: temperaturechart
                theme: ChartView.ChartThemeDark
                title: "Temperature"

                //anchors.fill: parent
                anchors.left: diffusionchart.right
                anchors.top: toptitle.bottom
                height: (root.height - toptitle.height)/2
                width: root.width/3
                antialiasing: true
                function updateETemp(x, y) {
                    temp_line.append(x,y);
                    if (x > temp_line.axisX.max) {
                        temp_line.axisX.max = x;
                    }
                    else if (x < temp_line.axisX.min) {
                        temp_line.axisX.min = x;
                    }
                    if (y > temp_line.axisY.max) {
                        temp_line.axisY.max = y;
                    }
                    else if (y < temp_line.axisY.min) {
                        temp_line.axisY.min = y;
                    }
                }
                function updateITemp(x, y) {
                    itemp_line.append(x,y);
                    if (x > itemp_line.axisX.max) {
                        itemp_line.axisX.max = x;
                    }
                    else if (x < itemp_line.axisX.min) {
                        itemp_line.axisX.min = x;
                    }
                    if (y > itemp_line.axisY.max) {
                        itemp_line.axisY.max = y;
                    }
                    else if (y < itemp_line.axisY.min) {
                        itemp_line.axisY.min = y;
                    }
                }

                LineSeries {
                    id: temp_line
                    axisX: ValueAxis {
                        min: 0
                        max: 0
                        tickCount: 5
                        titleText: "Step"
                    }

                    axisY: ValueAxis {
                        min: 0
                        max: 1e-7
                        titleText: "Temperature (K)"
                    }
                    name: "Emergent T"
                }
                LineSeries {
                    id: itemp_line
                    name: "Instantaneous T"
                }
            }
            ChartView {
                id: pressurechart
                theme: ChartView.ChartThemeDark
                title: "Pressure"

                //anchors.fill: parent
                anchors.left: parent.left
                anchors.top: energychart.bottom
                height: (root.height - toptitle.height)/2
                width: root.width/3
                antialiasing: true
                function updatePres(x, y) {
                    pres_line.append(x,y);
                    if (x > pres_line.axisX.max) {
                        pres_line.axisX.max = x;
                    }
                    else if (x < pres_line.axisX.min) {
                        pres_line.axisX.min = x;
                    }
                    if (y > pres_line.axisY.max) {
                        pres_line.axisY.max = y;
                    }
                    else if (y < pres_line.axisY.min) {
                        pres_line.axisY.min = y;
                    }
                }

                LineSeries {
                    id: pres_line
                    axisX: ValueAxis {
                        min: 0
                        max: 0
                        tickCount: 5
                        titleText: "Step"
                    }

                    axisY: ValueAxis {
                        min: 0
                        max: 1e-7
                        titleText: "I.G. pressure (atm)"
                    }
                    name: "Emergent Pressure"
                }

            }
            ChartView {
                id: velocitychart
                theme: ChartView.ChartThemeDark
                title: "Velocity"

                //anchors.fill: parent
                anchors.left: pressurechart.right
                anchors.top: diffusionchart.bottom
                height: (root.height - toptitle.height)/2
                width: root.width/3
                antialiasing: true
                function updateVel(x, y) {
                    vel_line.append(x,y);
                    if (x > vel_line.axisX.max) {
                        vel_line.axisX.max = x;
                    }
                    else if (x < vel_line.axisX.min) {
                        vel_line.axisX.min = x;
                    }
                    if (y > vel_line.axisY.max) {
                        vel_line.axisY.max = y;
                    }
                    else if (y < vel_line.axisY.min) {
                        vel_line.axisY.min = y;
                    }
                }

                LineSeries {
                    id: vel_line
                    axisX: ValueAxis {
                        min: 0
                        max: 0
                        tickCount: 5
                        titleText: "Step"
                    }

                    axisY: ValueAxis {
                        min: 0
                        max: 1e-7
                        titleText: "Velocity (m/s)"
                    }
                    name: "Average Velocity"
                }
            }


        }
        Page {
            Label {
                text: qsTr("Stuff on the 4 page")
            }
        }

        Page {
            Label {
                text: qsTr("more stuff, 5th")
            }
        }
    }

    footer: TabBar {
        id: tabBar
        currentIndex: swipeView.currentIndex
        TabButton {
            text: qsTr("Input setup")
        }
        TabButton {
            text: qsTr("MC")
        }
        TabButton {
            text: qsTr("MD")
        }
        TabButton {
            text: qsTr("hi")
        }
    }

}
