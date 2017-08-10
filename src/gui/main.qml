import QtQuick 2.7
import QtQuick.Controls 2.2
import QtQuick.Layouts 1.3
import QtQuick.Dialogs 1.2
import io.qt.examples.backend 1.0
import QtQuick.Controls.Styles 1.4
import QtQuick.Window 2.2
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
        source: "/home/khavernathy/mcmd/testzone/runlog.log"
        //linecount: 0 // this will change.
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
                        id: outputText
                        color: "white"
                        background: Rectangle { color: "black" }
                        width: 1000
                        height: 800
                        font.family: "Monospace"
                        //text: "Heloo world."
                        onTextChanged: {
                            console.log("the text changed yo")
                            //positionAt: bottom
                        }
                    }



            }

            Component.onCompleted: {
                //console.log( "WRITE"+ myFile.write("TEST"))
                outputText.text += runlogFile.read();
                flickableItem.contentY = flickableItem.contentHeight / 2 - height / 2
                flickableItem.contentX = flickableItem.contentWidth / 2 - width / 2

            }

            Button {
                id: reloadButton
                text: "reload"
                onClicked: {
                    outputText.text += runlogFile.read();
                    //outputText.text += "THE LINE COUNT IS: "+runlogFile.linecount;
                }
                anchors.right: parent.right
                anchors.bottom: parent.bottom

            }
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

        Page1 {


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
        TabButton {
            text: qsTr("testing ui")
        }
    }

}
