
#include "downloadManager.h"
#include <QFile>
#include <QFileInfo>
#include <QTextStream>


// this constructor get the url of a file and request the file to download
downloadManager::downloadManager(QUrl FileUrl, QObject *parent) : QObject(parent)
{
    QTextStream(stdout) << "downloadManager()" << endl;

    // Connect signal finished to the SLOT fileDownloaded()
    connect(&m_WebCtrl, SIGNAL (finished(QNetworkReply*)), this, SLOT ( fileDownloaded(QNetworkReply*) ) );

    QNetworkRequest request(FileUrl);
    m_WebCtrl.get(request);
}



downloadManager::~downloadManager() { }



// this slot it executed when the file is downloaded because it's connected befor in the constructor
void downloadManager::fileDownloaded(QNetworkReply* pReply) {

    QTextStream(stdout) << "fileDownloaded()" << endl;

    m_DownloadedData = pReply->readAll();

    QUrl url = pReply->url();
    if (pReply->error()) {
        fprintf(stderr, "Download of %s failed: %s\n", url.toEncoded().constData(), qPrintable(pReply->errorString()));
    }
    else {
        QString filename = saveFileName(url);
        if (saveToDisk(filename, pReply)) {
            printf("Download of %s succeeded (saved to %s)\n", url.toEncoded().constData(), qPrintable(filename));
        }
    }

    //emit a signal
    pReply->deleteLater();
    emit downloaded();

}



QByteArray downloadManager::downloadedData() const {

    QTextStream(stdout) << "downloadedData()" << endl;

    return m_DownloadedData;
}


// called from fileDownloaded()
QString downloadManager::saveFileName(const QUrl &url)
{
    QTextStream(stdout) << "saveFileName()" << endl;

    QString path = url.path();
    QString basename = QFileInfo(path).fileName();

    if (basename.isEmpty())
        basename = "download";

    if (QFile::exists(basename)) {
        // already exists, don't overwrite
        int i = 0;
        basename += '.';
        while (QFile::exists(basename + QString::number(i)))
            ++i;

        basename += QString::number(i);
    }

    return basename;
}


// called from fileDownloaded()
bool downloadManager::saveToDisk(const QString &filename, QIODevice *data)
{
    QTextStream(stdout) << "saveToDisk()" << endl;

    QFile file(filename);
    if (!file.open(QIODevice::WriteOnly)) {

        QTextStream(stdout) << "Could not open for writing" << endl;

        fprintf(stderr, "Could not open %s for writing: %s\n",
                qPrintable(filename),
                qPrintable(file.errorString()));
        return false;
    }

    file.write(data->readAll());

    QTextStream(stdout) << "file.write(data->readAll());" << endl;

    file.close();

    return true;
}
