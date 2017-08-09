import QtQuick 2.7
import QtQuick.Controls 2.2
import QtQuick.Layouts 1.3
import QtQuick.Dialogs 1.2
import io.qt.examples.backend 1.0
import FileIO 1.0



ApplicationWindow {
    id: root
    visible: true
    width: 1000
    height: 800
    title: qsTr("MCMD :: Monte Carlo & Molecular Dynamics")

    BackEnd {
        id: backend
        outputLineNumber: 0
    }
    FileIO {
        id: runlogFile
        source: "/Users/douglasfranz/mcmd/testzone/runlog.log"
        onError: console.log(msg);
    }
    SwipeView {
        id: swipeView
        anchors.fill: parent
        currentIndex: tabBar.currentIndex
        onCurrentIndexChanged: {
            //myText.text = myFile.read();

        }


        Page {

            ScrollView {
                id: view
                anchors.fill: parent
                TextArea {
                    width: 1000
                    height: 800
                    id: outputText
                    text: "Heloo world."
                    //anchors.centerIn: parent
                }
            }

            Component.onCompleted: {
                //console.log( "WRITE"+ myFile.write("TEST"))
                outputText.text = runlogFile.read();
            }

            Button {
                text: "reload"
                onClicked: outputText.text = runlogFile.read();

            }

            /*
            Label {
                text: qsTr("Page1 text")
            }
            TextField {
                text: backend.userName
                placeholderText: qsTr("User name")
                anchors.centerIn: parent
                onTextChanged: backend.userName = text
            }

            ListModel {
                id: model
                signal loadCompleted()

                Component.onCompleted: {
                    var xhr = new XMLHttpRequest;
                    xhr.open("GET", "/Users/douglasfranz/mcmd/examples/runlog.log");
                    xhr.onreadystatechange = function() {
                        if (xhr.readyState == XMLHttpRequest.DONE) {
                            var a = JSON.parse(xhr.responseText);
                            for (var b in a) {
                                var o = a[b];
                                model.append({label: o.label, value: o.value, description: o.description, image: o.image, selected: false});
                            }
                            model.loadCompleted()
                        }
                    }
                    xhr.send();
                }
            }

            TextField {
                height: 800
                width: 1000
                text: mcmddata.mcmdData
                placeholderText: qsTr("MCMD data will go here, i guess")+"<br />aldskjf"
            }
            */
        }

        Page {
            Label {
                text: qsTr("Second page")
                anchors.centerIn: parent
            }
            TextField {
                text: "Username: "+backend.userName+" ... and number: "+backend.outputLineNumber
            }
        }
        Page {
            Label {
                text: qsTr("Stuff on the 3rd page")
            }
        }
    }

    footer: TabBar {
        id: tabBar
        currentIndex: swipeView.currentIndex
        TabButton {
            text: qsTr("First")
        }
        TabButton {
            text: qsTr("Second")
        }
        TabButton {
            text: qsTr("3rd")
        }
    }
}
