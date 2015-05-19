import QtQuick 2.2
import QtQuick.Controls 1.1

ApplicationWindow {
    visible: true
    width: 640
    height: 480
    title: qsTr("Clicker")

    Column {
        id: leftColumn
        width: 300
        height: 100
        enabled: true

        Button {
            id: button1
            text: "pluss"
            width: 100
            height: 50
            anchors.left: parent.TopLeft
            property color buttonColor: "lightblue"
            property color onHoverColor: "gold"
            property color borderColor: "white"

        }

        MouseArea {
            id: mouseButton1
            anchors.fill: button1
            hoverEnabled: true
            onEntered: button1.border.color = onHoverColor
            onExited:  button1.border.color = borderColor
            property string plussbutton
            property int counter
            acceptedButtons: Qt.AllButtons

            onPressed: {
                if (mouse.button == Qt.LeftButton){
                    plussbutton = 'LeftMouseButton';
                    counter = counter + 1;
                }
                if (mouse.button == Qt.RightButton){
                    plussbutton = 'RightMouseButton'
                    counter = counter - 1;
                }
            }
        }


        Text {
            id: nClicks
            anchors.right: parent.right
            text: mouseButton1.counter
            font.pixelSize: 30
        }

    }

}
